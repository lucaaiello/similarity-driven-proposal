# Cold-chain postprocessing only; no compilation, fitting or sampling.
source("Real Data Analysis/Sampler/analysis_support.R")
apt_environment(); z <- long_load(); x <- z$inputs; p <- z$plan
out <- file.path(long_output, "summary")
apt_need(!file.exists(out), "Preserve existing summary")
paths <- unlist(lapply(1:4, function(start) file.path(long_output, long_methods,
  sprintf("ladder_%02d/checkpoints/iter_%06d", start, p$iterations))))
apt_need(all(dir.exists(paths)), "All12 final checkpoints must exist")
hashes <- tools::md5sum(file.path(paths, "state.rds"))
jp <- match("Prevotella", sub("^.*\\.", "", x$taxa))
ids <- match(c("3026", "4014"), rownames(x$X)); positive <- which(x$Y[, jp] > 0)
ix <- seq.int(p$burn + p$thin, p$iterations, p$thin)
post <- seq.int(p$burn + 1L, p$iterations)
records <- list(); cost <- list(); edges <- list(); pips <- list(); traces <- list(); locations <- list(); flows <- list()
for (i in seq_along(paths)) {
  s <- apt_read(paths[i]); method <- long_methods[(i-1L) %% 3L + 1L]; start <- (i-1L) %/% 3L + 1L
  expected <- p; expected$method <- method
  expected$seed_base <- expected$seed_base + match(method, long_methods) * 100L
  apt_need(identical(s$schema, "dm_apt50k_fresh_v1") && identical(s$plan, expected) &&
    identical(s$source_hashes, long_hashes()) && s$completed == p$iterations && s$ladder == start &&
    identical(s$input_hash, tools::md5sum("Real Data Analysis/Data/Prepared/analysis_plan.rds")) &&
    all(is.finite(s$ll_trace)) && all(is.finite(s$lp_trace)) &&
    length(s$cold_physical) == length(ix) && !any(vapply(s$cold_physical, is.null, TRUE)) &&
    all(colSums(s$swaps_attempted) == p$iterations) &&
    all(colSums(s$exchange_accepted) >= 0), "Checkpoint identity/completeness differs")
  pip <- Reduce("+", lapply(s$cold_physical, function(v) v$xi)) / length(ix)
  # Independent pure-R target reconstruction at the retained endpoints.
  for (d in c(1L, length(ix))) {
    ph <- s$cold_physical[[d]]; B <- matrix(0, ncol(x$X), ncol(x$Y))
    for (j in seq_len(ncol(x$Y))) B[ph$xi[, j] == 1L, j] <- ph$beta[[j]]
    target <- jb_target(ph$alpha, B, ph$xi, x$X, x$Y, rep(10, ncol(x$Y)), rep(10, ncol(x$Y)))
    apt_need(abs(target - s$ll_trace[ix[d], 1] - s$lp_trace[ix[d], 1]) < 1e-5,
      "Independent retained target reconstruction failed")
  }
  rec <- list(method = method, start = start, pips = pip, fitted = s$cold_fitted,
    model_size = s$scalar_trace[ix, 1, "model_size"], prev_size = s$scalar_trace[ix, 1, "prev_size"],
    prev_alpha = s$scalar_trace[ix, 1, "prev_alpha"],
    log_posterior = s$ll_trace[ix, 1] + s$lp_trace[ix, 1],
    subject3026 = s$cold_fitted[ids[1], ], subject4014 = s$cold_fitted[ids[2], ],
    positive_subject_mean = colMeans(s$cold_fitted[positive, , drop = FALSE]))
  records[[i]] <- rec
  cost[[i]] <- data.frame(method = method, start = long_names[start], retained = length(ix),
    seconds_including_startup = s$elapsed_seconds, allocated_core_hours = s$allocated_slot_seconds/3600,
    worker_cpu_hours_excluding_startup = sum(s$worker_cpu_seconds)/3600,
    cold_alpha_acceptance_full = s$counts[1, "alpha"]/p$iterations,
    cold_flip_acceptance_full = s$counts[1, "flip"]/(p$iterations*ncol(x$Y)),
    cold_local_accepted_full = s$counts[1, "local"], cold_joint_accepted_full = s$counts[1, "joint"],
    roundtrips_full = sum(s$flow$roundtrips), cold_swaps_post = s$swaps_accepted[2, 1],
    cold_swap_crossings4014_post = s$cold_swap_crossings[2, 2],
    cold_local_crossings4014_post = s$local_crossings[1, 2, 2])
  for (edge in 1:9) edges[[length(edges)+1L]] <- data.frame(method = method, start = long_names[start], edge = edge,
    beta_coldside_final = apt_beta(s$rho)[edge], beta_hotside_final = apt_beta(s$rho)[edge+1L],
    attempts_post = s$swaps_attempted[2, edge], accepted_post = s$swaps_accepted[2, edge],
    acceptance_post = s$swaps_accepted[2, edge]/s$swaps_attempted[2, edge],
    temperature_projections = s$rho_projection[edge])
  pg <- expand.grid(predictor = seq_len(ncol(x$X)), category = seq_len(ncol(x$Y)))
  pg$method <- method; pg$start <- long_names[start]; pg$pip <- as.vector(pip)
  pg$genus <- x$taxa[pg$category]; pg$predictor_code <- colnames(x$X)[pg$predictor]
  pips[[i]] <- pg
  traces[[i]] <- data.frame(method = method, start = long_names[start], iteration = ix,
    model_size = rec$model_size, prev_size = rec$prev_size, prev_alpha = rec$prev_alpha,
    log_posterior = rec$log_posterior, subject3026 = rec$subject3026, subject4014 = rec$subject4014)
  for (id in ids) locations[[length(locations)+1L]] <- data.frame(method = method, start = long_names[start],
    subject = rownames(x$X)[id], mean_percent = mean(s$cold_fitted[id, ])*100,
    fraction_below5percent = mean(s$cold_fitted[id, ] < .05))
  flows[[i]] <- data.frame(method = method, start = long_names[start], walker = 1:10,
    roundtrips = s$flow$roundtrips, distinct_temperature_slots = rowSums(s$visits > 0))
  rm(s); gc(FALSE)
}
bind <- function(v) { d <- do.call(rbind, v); rownames(d) <- NULL; d }
convergence <- list(); disagreement <- list()
for (method in long_methods) {
  rr <- records[vapply(records, function(r) r$method == method, TRUE)]
  for (field in c("model_size", "prev_size", "prev_alpha", "log_posterior", "subject3026", "subject4014", "positive_subject_mean")) {
    mat <- vapply(rr, function(r) r[[field]], numeric(length(ix)))
    convergence[[length(convergence)+1L]] <- data.frame(method = method, variable = field,
      rhat = jb_rhat(mat), range_of_chain_means = diff(range(colMeans(mat))))
  }
  pip <- simplify2array(lapply(rr, function(r) r$pips))
  span <- apply(pip, c(1, 2), function(v) diff(range(v)))
  fitted_span <- apply(vapply(rr, function(r) rowMeans(r$fitted), numeric(nrow(x$Y))), 1, function(v) diff(range(v)))
  disagreement[[length(disagreement)+1L]] <- data.frame(method = method, max_pip_spread = max(span),
    prev_pip_spreads_ge_half = sum(span[, jp] >= .5),
    max_positive_subject_spread_pp = max(fitted_span[positive])*100)
}
tables <- list(chain_metrics = bind(cost), exchange_edges = bind(edges), chain_pips = bind(pips),
  cold_traces = bind(traces), locations = bind(locations), walker_flow = bind(flows),
  convergence = bind(convergence), disagreement = bind(disagreement))
apt_need(identical(hashes, tools::md5sum(file.path(paths, "state.rds"))) && dir.create(out), "Inputs changed/output exists")
for (name in names(tables)) write.csv(tables[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE)
saveRDS(list(records = records, tables = tables, plan = p, checkpoint_hashes = hashes,
  source_hashes = long_hashes()), file.path(out, "summary.rds"))
report <- c("FRESH-START 50,000-ITERATION ADAPTIVE TEMPERING: ALL THREE SAMPLERS",
  "12 independent ladders: three samplers x four fresh null/5%/10%/15% starts.",
  "Ten replicas and ten one-thread workers per ladder. All replicas in each ladder share its fresh start.",
  "Burn25,000; no thinning; 25,000 cold draws per ladder. Hot draws do not enter posterior summaries.",
  "Every adjacent pair attempted once per iteration in random order, recalculating MH probabilities sequentially.",
  "Likelihood-only powers, unchanged priors/local kernels; MMV log-gap target0.234 and separate diversity adaptation online.",
  "All alpha/beta zero, g1 and alpha covariance identity initially; no old endpoint/tuning imported.",
  "Global alpha covariance is fixed apart from the existing temperature scaling; joint alpha/beta moves are enabled.",
  "Rhat/PIP disagreement must be assessed; longer chains or exchange acceptance alone do not establish convergence.",
  "Threshold crossings and walker roundtrips are diagnostic events, not independent posterior samples.",
  "Cross-method timing is descriptive: tasks may run on different hardware/load and uniform replicas can have different costs.",
  "", capture.output(print(tables$chain_metrics, row.names = FALSE)),
  "", capture.output(print(tables$convergence, row.names = FALSE)),
  "", capture.output(print(tables$disagreement, row.names = FALSE)))
writeLines(report, file.path(out, "report.txt")); writeLines("COMPLETE", file.path(out, "status.txt"))
cat(report, sep = "\n")
