# Local or cluster postprocessing only: no C++ compilation/optimizer/MCMC.
source("New_simulations/dm_adaptive_tempering/support.R")
source("New_simulations/dm_adaptive_tempering/reference.R")
apt_environment(); z <- apt_load(); x <- z$inputs; p <- z$plan
root <- "New_simulations/outputs/dm_apt10_v1"; out <- file.path(root, "summary")
apt_need(!file.exists(out), "Preserve existing APT summary")
paths <- c(file.path(root, sprintf("ladder_%02d/checkpoints/iter_%06d", 1:4, p$iterations)),
  file.path(root, sprintf("control_%02d/checkpoints/iter_%06d", 1:4, p$iterations)))
apt_need(all(dir.exists(paths)), "All four APT ladders and four matched controls must finish")
hashes <- tools::md5sum(file.path(paths, "state.rds"))
states <- lapply(paths, apt_read)
for (i in seq_along(states)) {
  s <- states[[i]]
  apt_need(s$completed == p$iterations && identical(s$plan, p) &&
    identical(s$source_hashes, apt_hashes()) && identical(s$input_hash, tools::md5sum(file.path(apt_root, "inputs/plan.rds"))) &&
    s$ladder == (i - 1L) %% 4L + 1L && identical(s$control, i > 4L), "APT task identity/provenance differs")
}
jp <- match("Prevotella", sub("^.*\\.", "", x$taxa)); positive <- which(x$Y[, jp] > 0)
ids <- match(c("3026", "4014"), rownames(x$X)); N <- p$iterations
ix <- seq.int(p$burn + p$thin, N, by = p$thin); post <- seq.int(p$burn + 1L, N)
records <- list(); cost <- list(); edges <- list(); flow <- list(); locations <- list(); levels <- list()
for (i in seq_along(states)) {
  s <- states[[i]]; arm <- if (s$control) "untempered_joint" else "adaptive_PT_joint"
  apt_need(length(s$cold_physical) == length(ix) && !any(vapply(s$cold_physical, is.null, TRUE)), "Cold draws incomplete")
  pip <- Reduce("+", lapply(s$cold_physical, function(v) v$xi)) / length(ix)
  for (d in c(1L, length(ix))) {
    ph <- s$cold_physical[[d]]; B <- matrix(0, ncol(x$X), ncol(x$Y))
    for (j in seq_len(ncol(x$Y))) B[ph$xi[, j] == 1L, j] <- ph$beta[[j]]
    target <- jb_target(ph$alpha, B, ph$xi, x$X, x$Y, rep(10, ncol(x$Y)), rep(10, ncol(x$Y)))
    apt_need(abs(target - (s$ll_trace[ix[d], 1] + s$lp_trace[ix[d], 1])) < 1e-5,
      "Independent cold endpoint target reconstruction failed")
  }
  records[[i]] <- list(arm = arm, ladder = s$ladder, pips = pip,
    alpha = vapply(s$cold_physical, function(v) v$alpha[jp], 0),
    prev_size = vapply(s$cold_physical, function(v) sum(v$xi[, jp]), 0),
    size = vapply(s$cold_physical, function(v) sum(v$xi), 0), fitted = s$cold_fitted,
    posterior = s$ll_trace[ix, 1] + s$lp_trace[ix, 1],
    positive_mean = colMeans(s$cold_fitted[positive, , drop = FALSE]))
  changes <- sum(vapply(2:length(s$cold_physical), function(d)
    any(s$cold_physical[[d]]$xi[, jp] != s$cold_physical[[d - 1L]]$xi[, jp]), TRUE))
  cost[[i]] <- data.frame(arm = arm, ladder = s$ladder, replicas = length(s$slots),
    elapsed_seconds_including_worker_compilation = s$elapsed_seconds,
    allocated_core_hours = s$allocated_slot_seconds / 3600,
    propagation_worker_cpu_hours_excluding_compilation = sum(s$worker_cpu_seconds) / 3600,
    observed_prev_configuration_changes_thinned = changes,
    cold_local_crossings4014_post = s$local_crossings[1, 2, 2],
    cold_swap_crossings4014_post = s$cold_swap_crossings[2, 2],
    roundtrips_full_run = if (s$control) NA_integer_ else sum(s$flow$roundtrips),
    roundtrips_completed_after_burn = if (s$control) NA_integer_ else tail(s$roundtrip_trace, 1) - s$roundtrip_trace[p$burn])
  if (!s$control) {
    L <- length(s$slots); probabilities <- colMeans(tail(s$hypothetical_swap, 500))
    for (edge in seq_len(L - 1L)) edges[[length(edges) + 1L]] <- data.frame(ladder = s$ladder,
      edge = edge, beta_coldside_final = apt_beta(s$rho)[edge], beta_hotside_final = apt_beta(s$rho)[edge + 1L],
      attempts_post = s$swaps_attempted[2, edge], accepted_post = s$swaps_accepted[2, edge],
      acceptance_post = if (s$swaps_attempted[2, edge]) s$swaps_accepted[2, edge] / s$swaps_attempted[2, edge] else NA_real_,
      mean_hypothetical_acceptance_last500 = probabilities[edge], projected_updates = s$rho_projection[edge])
    flow[[length(flow) + 1L]] <- data.frame(ladder = s$ladder, walker = seq_len(L),
      initial_origin = s$seed_origin, roundtrips = s$flow$roundtrips,
      distinct_temperature_slots_visited = rowSums(s$visits > 0))
    for (level in seq_len(L)) levels[[length(levels) + 1L]] <- data.frame(ladder = s$ladder,
      level = level, inverse_temperature_final = apt_beta(s$rho)[level],
      local_crossings3026_post = s$local_crossings[level, 1, 2],
      local_crossings4014_post = s$local_crossings[level, 2, 2],
      fraction4014_below5_post = mean(s$scalar_trace[post, level, "subject4014"] < .05),
      fraction4014_above20_post = mean(s$scalar_trace[post, level, "subject4014"] > .2),
      ll_q10_last500 = quantile(tail(s$ll_trace[, level], 500), .1, names = FALSE),
      ll_q50_last500 = median(tail(s$ll_trace[, level], 500)),
      ll_q90_last500 = quantile(tail(s$ll_trace[, level], 500), .9, names = FALSE),
      g_mean_final = mean(s$slots[[level]]$tuning$g))
  }
  for (id in ids) locations[[length(locations) + 1L]] <- data.frame(arm = arm, ladder = s$ladder,
    subject = rownames(x$X)[id], mean_percent = mean(s$cold_fitted[id, ]) * 100,
    min_percent = min(s$cold_fitted[id, ]) * 100, max_percent = max(s$cold_fitted[id, ]) * 100,
    proportion_below5percent = mean(s$cold_fitted[id, ] < .05),
    proportion_above20percent = mean(s$cold_fitted[id, ] > .2))
}
bind <- function(v) { d <- do.call(rbind, v); rownames(d) <- NULL; d }
cost <- bind(cost); edges <- bind(edges); flow <- bind(flow); locations <- bind(locations); levels <- bind(levels)
convergence <- list(); disagreement <- list()
for (arm in unique(cost$arm)) {
  rr <- records[vapply(records, function(r) r$arm == arm, TRUE)]
  fields <- list(prev_alpha = vapply(rr, function(r) r$alpha, numeric(length(ix))),
    prev_size = vapply(rr, function(r) r$prev_size, numeric(length(ix))),
    model_size = vapply(rr, function(r) r$size, numeric(length(ix))),
    posterior = vapply(rr, function(r) r$posterior, numeric(length(ix))),
    subject3026 = vapply(rr, function(r) r$fitted[ids[1], ], numeric(length(ix))),
    subject4014 = vapply(rr, function(r) r$fitted[ids[2], ], numeric(length(ix))),
    positive_subject_mean = vapply(rr, function(r) r$positive_mean, numeric(length(ix))))
  for (name in names(fields)) convergence[[length(convergence) + 1L]] <- data.frame(arm = arm,
    variable = name, rhat = jb_rhat(fields[[name]]),
    constant_chains = sum(apply(fields[[name]], 2, function(v) length(unique(v)) == 1L)),
    chain_mean_range = diff(range(colMeans(fields[[name]]))))
  ps <- simplify2array(lapply(rr, function(r) r$pips)); spans <- apply(ps, c(1, 2), function(v) diff(range(v)))
  fitted_means <- vapply(rr, function(r) rowMeans(r$fitted), numeric(nrow(x$Y)))
  fspan <- apply(fitted_means, 1, function(v) diff(range(v))) * 100
  disagreement[[length(disagreement) + 1L]] <- data.frame(arm = arm, max_pip_spread = max(spans),
    prev_spreads_ge_half = sum(spans[, jp] >= .5), positive_subjects_spread_gt5pp = sum(fspan[positive] > 5),
    positive_subject_median_spread_pp = median(fspan[positive]), positive_subject_max_spread_pp = max(fspan[positive]))
}
convergence <- bind(convergence); disagreement <- bind(disagreement)
apt_need(identical(hashes, tools::md5sum(file.path(paths, "state.rds"))) && dir.create(out), "Inputs changed or summary cannot be created")
saveRDS(list(records = records, costs = cost, edges = edges, flow = flow, locations = locations, levels = levels,
  convergence = convergence, disagreement = disagreement, plan = p, source_hashes = apt_hashes(),
  checkpoint_hashes = hashes, temperature_histories = lapply(states[1:4], function(s) s$beta_trace),
  walker_histories = lapply(states[1:4], function(s) s$walker_trace)), file.path(out, "summary.rds"))
fmt <- function(v) capture.output(print(v, row.names = FALSE, digits = 5))
report <- c("MMV ADAPTIVE LIKELIHOOD-TEMPERING FEASIBILITY TEST",
  "4 independent ten-replica ladders plus4 matched untempered joint-update controls.",
  "5000 sweeps per replica;first2000 excluded;thin5;online temperature adaptation targets0.234.",
  "MMV log-gap adaptation is applied to the likelihood-only path, not the complete published fixed-dimension RW algorithm.",
  "The cold target is the unchanged posterior. Hot-replica draws are not used for scientific posterior summaries.",
  "The existing effective-neighbor/diversity adaptation remains online and separate from the swap target.",
  "Initial hot replicas use four KNOWN saved-origin configurations. Imported seeded regions do not prove discovery of new modes.",
  "Count both actual hot LOCAL transitions and cold exchanges;roundtrips trace physical walkers, not independent posterior samples.",
  "All reported crossings use fitted-fraction5% as a descriptive diagnostic, not a selection threshold.",
  "Four independent ladders supply four cold chains;the10 temperatures are not independent posterior replications.",
  "Paper convergence/optimality theorems do not automatically transfer to this trans-dimensional likelihood-path implementation.",
  "No automatic claim that10 replicas or0.234 exchange acceptance is sufficient is made.",
  "", "COST AND COLD-CHAIN TRANSITIONS", fmt(cost), "", "ADJACENT EDGES", fmt(edges),
  "", "WALKER FLOW (full run, includes temperature adaptation)", fmt(flow),
  "", "LOCAL EXPLORATION AT EACH TEMPERATURE (post-burn)", fmt(levels),
  "", "COLD FITTED FRACTIONS AFTER BURN", fmt(locations),
  "", "COLD-CHAIN FOLDED RANK RHAT", fmt(convergence),
  "", "COLD-CHAIN POSTERIOR DISAGREEMENT", fmt(disagreement),
  "", "SUCCESS REQUIRES:hot local crossings,ladder communication,and reduced cold posterior disagreement per computational budget.",
  "Acceptance near0.234 alone is not success. This bounded feasibility test does not certify converged biological PIPs.")
writeLines(report, file.path(out, "report.txt")); writeLines("COMPLETE", file.path(out, "status.txt")); cat(report, sep = "\n")
