# parse/compile/fixed modes never run MCMC or an optimizer.
# full mode adds bounded SYNTHETIC sampler/parallel/restart checks only.
source("New_simulations/dm_adaptive_tempering/support.R")
source("New_simulations/dm_adaptive_tempering/math.R")
source("New_simulations/dm_adaptive_tempering/reference.R")
source("New_simulations/dm_adaptive_tempering/runner.R")
apt_environment()
o <- apt_options(c("mode", "output")); mode <- apt_get(o, "mode", "full")
apt_need(mode %in% c("parse", "compile", "fixed", "full"), "Invalid audit mode")
for (p in apt_sources()[grepl("[.]R$", apt_sources())]) invisible(parse(p))
if (mode == "parse") { cat("APT PARSE PASSED; no compiler, optimizer or MCMC\n"); quit(save = "no") }
out <- apt_get(o, "output", file.path("New_simulations/outputs", paste0("dm_apt_audit_",
  format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())))
apt_need(!file.exists(out) && dir.create(out, recursive = TRUE), "Preserve existing audit output")
writeLines("INCOMPLETE", file.path(out, "status.txt")); sources <- apt_hashes()
checks <- list()
check <- function(name, ok, detail = "") {
  checks[[length(checks) + 1L]] <<- list(check = name, passed = isTRUE(ok), detail = detail)
  if (!isTRUE(ok)) {
    saveRDS(checks, file.path(out, "checks.rds")); writeLines("FAILED", file.path(out, "status.txt"))
    stop("APT audit failed: ", name, " ", detail, call. = FALSE)
  }
}
near <- function(a, b, tol = 1e-9) isTRUE(all.equal(unname(a), unname(b), tolerance = tol, check.attributes = FALSE))
rho <- rep(log(log(32) / 9), 9)
b <- apt_beta(rho)
check("MMV_initial_geometric_ladder", near(b, exp(-seq(0, log(32), length.out = 10))))
check("MMV_cold_exactly_one_ordered_positive", b[1] == 1 && all(diff(b) < 0) && min(b) > 0)
u <- apt_rho_update(rho, rep(0, 10), 3)
check("MMV_high_acceptance_widens_every_gap", all(u$rho > rho))
ll <- seq(0, -1e9, length.out = 10)
u <- apt_rho_update(rho, ll, 3)
check("MMV_low_acceptance_narrows_every_gap", all(u$rho < rho))
check("MMV_all_pairs_not_binary_accepted_flags", near(u$rho, rho + 4^(-.6) * (apt_swap_prob(ll, b) - .234)))
check("MMV_projection_enforced", all(apt_rho_update(rep(log(2), 9), rep(0, 10), 1)$rho <= log(2)))
check("MMV_probability_logstable", all(apt_swap_prob(seq(0, -1e300, length.out = 10), b) == 0))
check("MMV_diminishing_gain", apt_rho_update(rho, ll, 100)$gain < apt_rho_update(rho, ll, 10)$gain)

# Exact finite-state stationarity checks, no random simulation. Unequal priors
# are deliberate: likelihood tempering must not flatten the model prior.
lll <- c(-7, -4, 1, .5); prior <- c(.1, .5, .3, .1); beta <- c(1, .4, .1)
grid <- as.matrix(expand.grid(1:4, 1:4, 1:4)); K <- nrow(grid)
joint <- exp(rowSums(matrix(lll[grid], K, 3) * rep(beta, each = K)) + rowSums(matrix(log(prior)[grid], K, 3)))
joint <- joint / sum(joint); key <- apply(grid, 1, paste, collapse = "/")
Q <- matrix(0, K, K)
for (i in seq_len(K)) for (edge in 1:2) {
  next_state <- grid[i, ]; next_state[c(edge, edge + 1L)] <- next_state[c(edge + 1L, edge)]
  j <- match(paste(next_state, collapse = "/"), key)
  a <- apt_swap_prob(lll[grid[i, ]], beta)[edge]
  Q[i, j] <- Q[i, j] + a / 2; Q[i, i] <- Q[i, i] + (1 - a) / 2
}
check("finite_exchange_rows_sum_one", max(abs(rowSums(Q) - 1)) < 1e-14)
check("finite_exchange_detailed_balance", max(abs(joint * Q - t(joint * Q))) < 1e-14)
check("finite_exchange_stationary", max(abs(drop(joint %*% Q) - joint)) < 1e-14)
cold <- vapply(1:4, function(v) sum(joint[grid[, 1] == v]), 0)
target <- prior * exp(lll); target <- target / sum(target)
check("finite_cold_marginal_original_posterior", max(abs(cold - target)) < 1e-14)
for (t in beta) {
  w <- prior * exp(t * lll); w <- w / sum(w); M <- matrix(0, 4, 4)
  for (i in 1:4) for (j in 1:4) {
    a <- min(1, w[j] / w[i]); M[i, j] <- M[i, j] + a / 4
    M[i, i] <- M[i, i] + (1 - a) / 4
  }
  check(paste("finite_local_tempered_stationary", t), max(abs(drop(w %*% M) - w)) < 1e-14)
}
new <- apt_compile(FALSE); baseline <- apt_compile(TRUE)
if (mode == "compile") {
  writeLines("COMPILE_ONLY", file.path(out, "status.txt")); cat("APT COMPILE PASSED; no MCMC or optimizer\n"); quit(save = "no")
}
set.seed(2026102801)
X <- unclass(scale(matrix(rnorm(48), 12, 4))); rownames(X) <- c("3026", "4014", as.character(3:12))
Y <- t(vapply(1:12, function(i) as.integer(rmultinom(1, 35, c(.2, .3, .5))), integer(3)))
for (t in c(1, .4, .01)) for (k in c(0L, 1L, 3L)) for (shift in c(0, 40, 1000)) {
  alpha <- c(-.3, .2, .1) + shift; xi <- matrix(0L, 4, 3); B <- matrix(0, 4, 3)
  if (k) { xi[1:k, 2] <- 1L; B[1:k, 2] <- seq_len(k) * .07 }
  selected <- which(xi[, 2] == 1L); design <- X[, selected, drop = FALSE]; beta <- B[selected, 2]
  candidate <- c(alpha[2] + .05, beta + .02 * seq_along(beta))
  ca <- alpha; ca[2] <- candidate[1]; cb <- B; cb[selected, 2] <- candidate[-1]
  s2 <- c(5, 7, 9); r2 <- c(10, 11, 12)
  ref_ll <- jb_likelihood(Y, sweep(X %*% B, 2, alpha, "+"))
  next_ll <- jb_likelihood(Y, sweep(X %*% cb, 2, ca, "+"))
  dp <- -.5 * (candidate[1]^2 - alpha[2]^2) / s2[2] - .5 * (sum(candidate[-1]^2) - sum(beta^2)) / r2[2]
  raw <- jb_raw_precision(design, Y, alpha, X %*% B, 2, beta, s2[2], r2[2])
  rawr <- jb_raw_precision(design, Y, ca, X %*% cb, 2, candidate[-1], s2[2], r2[2])
  priorQ <- diag(c(1 / s2[2], rep(1 / r2[2], k)), k + 1L)
  raw <- t * (raw - priorQ) + priorQ; rawr <- t * (rawr - priorQ) + priorQ
  qf <- jb_project(raw) * (k + 1) / 4; qr <- jb_project(rawr) * (k + 1) / 4
  q <- new$audit_dm_joint_block_terms(design, Y, alpha, X %*% B, 2L, beta, candidate, s2[2], r2[2], 2, t)
  label <- paste(t, k, shift, sep = "/")
  check(paste(label, "joint_forward_actual_prior_precision"), near(q$forward_precision, qf))
  check(paste(label, "joint_reverse_actual_prior_precision"), near(q$reverse_precision, qr))
  expected <- t * (next_ll - ref_ll) + dp + jb_logq(c(alpha[2], beta), candidate, qr) - jb_logq(candidate, c(alpha[2], beta), qf)
  check(paste(label, "joint_tempered_MH"), abs(q$ratio - expected) < 5e-7)
  reverse <- new$audit_dm_joint_block_terms(design, Y, ca, X %*% cb, 2L, candidate[-1], c(alpha[2], beta), s2[2], r2[2], 2, t)
  check(paste(label, "joint_antisymmetry"), abs(q$ratio + reverse$ratio) < 1e-7)
  check(paste(label, "prior_and_q_are_not_tempered"),
    abs(new$audit_dm_pt_ratio(t, next_ll - ref_ll, dp, 3.7) - (t * (next_ll - ref_ll) + dp + 3.7)) < 1e-10)
  if (k) for (fixed in c(FALSE, TRUE)) {
    cn <- beta + .025; nb <- B; nb[selected, 2] <- cn
    w <- new$audit_dm_within_model_terms(design, Y, alpha, X %*% B, 2L, beta, cn, r2[2], .25, fixed, t)
    dll <- jb_likelihood(Y, sweep(X %*% nb, 2, alpha, "+")) - ref_ll
    expected <- t * dll - .5 * (sum(cn^2) - sum(beta^2)) / r2[2] + if (fixed) 0 else w$log_reverse - w$log_forward
    check(paste(label, "beta_within_tempered", fixed), abs(w$ratio - expected) < 5e-7)
  }
}

# Full-state exchanges: dimension-varying beta and indicators travel together;
# RNG and temperature-slot tuning are intentionally not exchanged.
fake <- list(labels = 1:2, slots = list(list(physical = list(beta = list(1), xi = 1), ll = 2,
  lp = 3, score = list(v = 4), rng = 11, tuning = 12), list(physical = list(beta = list(1:2), xi = 2),
  ll = 5, lp = 6, score = list(v = 7), rng = 21, tuning = 22)))
swapped <- apt_swap_slots(fake, 1L)
check("exchange_complete_physical_state", identical(swapped$slots[[1]]$physical, fake$slots[[2]]$physical))
check("exchange_keeps_temperature_tuning_and_rng", swapped$slots[[1]]$rng == 11 && swapped$slots[[1]]$tuning == 12)
check("exchange_reverse_involution", identical(apt_swap_slots(swapped, 1L), fake))
flow <- list(phase = c(1L, 0L), departed = c(0L, 0L), roundtrips = c(0L, 0L), durations = integer())
flow <- apt_flow_update(flow, 2:1, 10L); flow <- apt_flow_update(flow, 1:2, 20L)
check("walker_cold_hot_cold_roundtrip", flow$roundtrips[1] == 1L && flow$durations[1] == 20L)

if (mode == "full") {
  z <- apt_load(); plan <- z$plan
  x <- list(X = X, Y = Y, neighbors = lapply(1:4, function(j) as.integer(setdiff(1:4, j))),
    taxa = c("Bacteroides", "Prevotella", "Control"), origins = list())
  for (chain in 1:4) x$origins[[sprintf("online_similarity_local_chain_%02d", chain)]] <- list(
    start = list(alpha = rep(0, 3), xi = matrix(0L, 4, 3), beta = rep(list(numeric()), 3),
      g = rep(1, 3), Sigma_alpha = diag(3) * .01))
  assign(".apt_x", x, .GlobalEnv); assign(".apt_cpp", new, .GlobalEnv)
  initial <- apt_worker_initial_scores("online_similarity_local")
  frozen <- plan; frozen$adapt_g <- FALSE; frozen$replicas <- 1L; frozen$iterations <- 30L
  frozen$burn <- 0L; frozen$thin <- 1L
  single <- apt_new_state(x, frozen, 1L, TRUE, initial)
  args <- apt_args(x, single$slots[[1]]$physical, single$slots[[1]]$tuning, 1, 30L)
  old_args <- args[setdiff(names(args), c("inverse_temperature", "quiet"))]
  assign(".Random.seed", single$slots[[1]]$rng, .GlobalEnv)
  old <- do.call(baseline$DM_joint_block_VS_RJMCMC, old_args); old_rng <- .Random.seed
  serial <- function(requests) lapply(requests, apt_worker_step)
  for (i in 1:30) single <- apt_iteration(single, x, serial)
  check("cold_30sweeps_RNG_original", identical(single$slots[[1]]$rng, old_rng))
  check("cold_30sweeps_alpha_original", near(single$slots[[1]]$physical$alpha, tail(old$alpha, 1)))
  check("cold_30sweeps_indicator_original", identical(single$slots[[1]]$physical$xi, old$xi[, , 30]))
  check("cold_30sweeps_beta_original", near(single$slots[[1]]$physical$beta, old$beta[[30]]))
  for (t in c(1, .4, .01)) {
    slot <- single$slots[[1]]
    for (i in 1:10) {
      v <- apt_worker_step(list(physical = slot$physical, tuning = slot$tuning, rng = slot$rng,
        beta = t, method = "online_similarity_local"))
      pp <- v$physical; BB <- matrix(0, 4, 3)
      for (j in 1:3) BB[pp$xi[, j] == 1L, j] <- pp$beta[[j]]
      ref <- jb_target(pp$alpha, BB, pp$xi, X, Y, rep(10, 3), rep(10, 3))
      check(paste("hot_independent_retained_target", t, i), abs(v$ll + v$lp - ref) < 1e-7)
      slot$physical <- v$physical; slot$rng <- v$rng
    }
  }
  # External online diversity update against the unchanged C++ rule.
  adaptive <- old_args; adaptive$adapt_g <- TRUE; adaptive$adaptation_start <- 1L
  adaptive$adaptation_end <- 30L; adaptive$window_size <- 5L
  assign(".Random.seed", apt_new_state(x, frozen, 1L, TRUE, initial)$slots[[1]]$rng, .GlobalEnv)
  old_adapt <- do.call(baseline$DM_joint_block_VS_RJMCMC, c(adaptive[names(adaptive) != "pilot_start"],
    list(pilot_start = NULL)))
  # Compare the formula directly using recorded signals, independent of the
  # different cold-start intercept covariance used above.
  tuning <- list(g = rep(1, 3), effective_sum = rep(0, 3), effective_count = 0L, g_epoch = 0L)
  gp <- plan; gp$g_window <- 5L; gp$g_target <- 4
  for (i in 1:30) {
    tuning <- apt_g_update(tuning, as.numeric(old_adapt$proposal_effective_neighbors_trace[i, ]), gp)
    check(paste("online_g_Cpp_formula", i), near(tuning$g, as.numeric(old_adapt$g_trace[i, ])))
  }
  short <- plan; short$replicas <- 3L; short$iterations <- 24L; short$burn <- 8L
  short$thin <- 2L; short$g_window <- 5L; short$g_target <- 4; short$adapt_until <- 24L
  start <- apt_new_state(x, short, 2L, FALSE, initial)
  full <- start
  for (i in 1:24) full <- apt_iteration(full, x, serial)
  part <- start
  for (i in 1:9) part <- apt_iteration(part, x, serial)
  part <- unserialize(serialize(part, NULL))
  for (i in 10:24) part <- apt_iteration(part, x, serial)
  fields <- setdiff(names(full), c("native_seconds", "worker_cpu_seconds"))
  check("serialized_resume_all_physical_and_adaptive_state", near(full[fields], part[fields], 1e-12))
  check("serialized_resume_exact_walker_and_rng", identical(full$labels, part$labels) &&
    identical(full$master_rng, part$master_rng) && all(vapply(1:3, function(k) identical(full$slots[[k]]$rng, part$slots[[k]]$rng), TRUE)))
  requests <- lapply(1:3, function(k) list(physical = start$slots[[k]]$physical, tuning = start$slots[[k]]$tuning,
    rng = start$slots[[k]]$rng, beta = apt_beta(start$rho)[k], method = plan$method))
  expected <- serial(requests)
  cl <- apt_pool(x, 2L)
  actual <- tryCatch(parallel::parLapply(cl, requests, apt_worker_step), finally = parallel::stopCluster(cl))
  for (i in 1:3) {
    check(paste("PSOCK_exact_rng_and_indicators", i), identical(expected[[i]]$rng, actual[[i]]$rng) &&
      identical(expected[[i]]$physical$xi, actual[[i]]$physical$xi))
    f <- setdiff(names(expected[[i]]), c("local_seconds", "cpu_seconds"))
    check(paste("PSOCK_target_and_state", i), near(expected[[i]][f], actual[[i]][f]))
  }
}
apt_need(identical(sources, apt_hashes()), "APT sources changed during audit")
saveRDS(checks, file.path(out, "checks.rds"))
if (mode == "full") saveRDS(list(passed = TRUE, source_hashes = sources,
  input_hash = tools::md5sum(file.path(apt_root, "inputs/plan.rds")), runtime = apt_runtime(), checks = length(checks)),
  file.path(out, "audit.rds"))
report <- c(paste("ADAPTIVE TEMPERING", toupper(mode), "AUDIT PASSED"), paste("Checks:", length(checks)),
  "MMV log-gap adaptation, likelihood-only powers, untempered priors and reverse Gaussian geometry checked.",
  "Exact finite-state swap invariance/cold marginal and full-state/tuning separation checked without MCMC.",
  if (mode == "full") "Bounded SYNTHETIC MCMC, original cold kernel, parallel RNG and serialized restart checked. No real-data pilot was run."
  else "No sampler or optimizer executed. Fixed-only success does not permit a real-data pilot.")
writeLines(report, file.path(out, "report.txt")); cat(report, sep = "\n")
writeLines(if (mode == "full") "COMPLETE" else "FIXED_ONLY", file.path(out, "status.txt"))
