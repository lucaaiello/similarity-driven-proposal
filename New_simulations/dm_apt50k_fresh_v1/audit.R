# Pure finite-state checks plus bounded SYNTHETIC MCMC only; never real production.
source("New_simulations/dm_apt50k_fresh_v1/support.R")
source("New_simulations/dm_apt50k_fresh_v1/runner.R")
apt_environment(); o <- apt_options(c("output", "workers"))
workers <- as.integer(apt_get(o, "workers", "2"))
apt_need(workers %in% c(2L, 10L), "Audit workers must be2 or10")
out <- apt_get(o, "output", file.path("New_simulations/outputs", paste0("dm_apt50k_audit_",
  format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())))
apt_need(!file.exists(out) && dir.create(out, recursive = TRUE), "Preserve existing audit output")
writeLines("INCOMPLETE", file.path(out, "status.txt"))
hash <- long_hashes(); protected <- aptopt_protected(); z <- long_load(); checks <- list()
check <- function(name, ok) {
  checks[[length(checks) + 1L]] <<- list(check = name, passed = isTRUE(ok))
  if (!isTRUE(ok)) {
    saveRDS(checks, file.path(out, "checks.rds")); writeLines("FAILED", file.path(out, "status.txt"))
    stop("Fresh50K audit failed: ", name, call. = FALSE)
  }
}
# Independent exact finite-state composition: unequal priors, likelihood powers.
ll <- c(-3, .5, 2); prior <- c(.2, .5, .3); beta <- c(1, .4, .1)
grid <- as.matrix(expand.grid(1:3, 1:3, 1:3)); K <- nrow(grid)
key <- apply(grid, 1, paste, collapse = "/")
pi <- exp(rowSums(matrix(ll[grid], K, 3) * rep(beta, each = K)) +
  rowSums(matrix(log(prior)[grid], K, 3))); pi <- pi / sum(pi)
edges <- lapply(1:2, function(edge) {
  M <- matrix(0, K, K)
  for (i in 1:K) {
    target <- grid[i, ]; target[c(edge, edge+1L)] <- target[c(edge+1L, edge)]
    j <- match(paste(target, collapse = "/"), key)
    a <- min(1, exp((beta[edge] - beta[edge+1L]) * (ll[grid[i, edge+1L]] - ll[grid[i, edge]])))
    M[i, j] <- M[i, j] + a; M[i, i] <- M[i, i] + 1-a
  }
  M
})
Q <- (edges[[1]] %*% edges[[2]] + edges[[2]] %*% edges[[1]]) / 2
check("random_permutation_rows", max(abs(rowSums(Q)-1)) < 1e-14)
check("random_permutation_joint_stationary", max(abs(drop(pi %*% Q)-pi)) < 1e-14)
check("random_permutation_reversible_mixture", max(abs(pi * Q - t(pi * Q))) < 1e-14)
cold <- vapply(1:3, function(k) sum(pi[grid[, 1] == k]), 0)
target <- prior * exp(ll); target <- target / sum(target)
check("cold_marginal_preserved", max(abs(cold-target)) < 1e-14)
for (method in long_methods) for (start in 1:4) {
  v <- z$inputs$origins[[sprintf("%s_chain_%02d", method, start)]]
  st <- v$start
  check(paste("fresh_inputs", method, start), v$origin_iteration == 0L && all(st$alpha == 0) &&
    all(unlist(st$beta) == 0) && all(st$g == 1) && identical(st$Sigma_alpha, diag(ncol(z$inputs$Y))) &&
    sum(st$xi) == round(prod(dim(st$xi)) * c(0, .05, .10, .15)[start]))
}
cat("BOUNDed SYNTHETIC audit: local sweeps, 10-replica exchanges, adaptation, parallel restart.\n")
x <- long_fresh(aptopt_synthetic_fixture())
aptopt_initialize(x, "both")
p <- z$plan; p$iterations <- 24L; p$burn <- 12L; p$thin <- 1L
p$adapt_until <- 24L; p$g_window <- 5L; p$g_target <- 5
run <- function(method, engine, cl = NULL, restart = FALSE) {
  p$method <- method
  initial <- long_seed_scores(x, method, 2L)
  s <- long_new(x, p, 2L, initial)
  prop <- function(requests) {
    requests <- lapply(requests, function(r) { r$engine <- engine; r })
    if (is.null(cl)) lapply(requests, aptopt_step) else parallel::parLapply(cl, requests, aptopt_step)
  }
  for (n in 1:p$iterations) {
    s <- long_iteration(s, x, prop)
    if (restart && n == 11L) {
      saveRDS(s, file.path(out, paste0(method, "_restart.rds")))
      s <- readRDS(file.path(out, paste0(method, "_restart.rds")))
      if (is.null(cl)) aptopt_reset_contexts() else parallel::clusterCall(cl, aptopt_reset_contexts)
    }
  }
  s
}
cl <- aptopt_pool(x, workers, "persistent")
tryCatch({
  for (method in long_methods) {
    a <- run(method, "reference")
    b <- run(method, "persistent", cl, TRUE)
    check(paste("full_ladder_parallel_engine_restart", method), aptopt_ladder_comparison(a, b)$passed &&
      identical(a$exchange_order, b$exchange_order) && identical(a$exchange_accepted, b$exchange_accepted))
    check(paste("each_edge_once_every_sweep", method), all(apply(a$exchange_order, 1,
      function(v) identical(sort(v), 1:9))) && sum(a$swaps_attempted) == 24L * 9L &&
      all(colSums(a$swaps_attempted) == 24L))
    check(paste("accepted_swaps_and_online_g", method), sum(a$swaps_accepted) > 0 &&
      all(vapply(a$slots, function(v) if (method == "uniform_flip") all(v$tuning$g == 1) else v$tuning$g_epoch > 0, TRUE)))
    check(paste("all_replicas_own_start", method), all(a$seed_origin == 2L))
    check(paste("unthinned_retained_draws", method), length(a$cold_physical) == 12L &&
      !any(vapply(a$cold_physical, is.null, TRUE)) &&
      identical(as.numeric(a$scalar_trace[13:24, 1, "model_size"]),
        vapply(a$cold_physical, function(ph) as.numeric(sum(ph$xi)), 0)))
  }
}, finally = parallel::stopCluster(cl))
# A forced accepted walk through the actual production exchange code checks
# sequential probability recomputation, full-state exchange and slot RNG/tuning.
fake <- long_new(x, p, 2L, long_seed_scores(x, p$method, 2L))
for (k in 1:10) { fake$slots[[k]]$ll <- k * .2; fake$slots[[k]]$lp <- k; }
set.seed(8723); saved <- .Random.seed; expected <- fake
order <- sample.int(9L)
for (pos in seq_along(order)) {
  e <- order[pos]; beta <- apt_beta(fake$rho)
  prob <- min(1, exp((beta[e]-beta[e+1L]) *
    (expected$slots[[e+1L]]$ll - expected$slots[[e]]$ll)))
  if (runif(1) < prob) expected <- apt_swap_slots(expected, e)
}
expected_rng <- .Random.seed
assign(".Random.seed", saved, .GlobalEnv)
got <- long_exchange(fake, x, 1L, apt_beta(fake$rho))
check("actual_exchange_sequence_reference", identical(got$slots, expected$slots) &&
  identical(got$labels, expected$labels) && identical(.Random.seed, expected_rng))
check("protected_originals_unchanged", identical(protected, aptopt_protected()))
check("new_sources_unchanged", identical(hash, long_hashes()))
saveRDS(checks, file.path(out, "checks.rds"))
saveRDS(list(passed = TRUE, source_hashes = hash, runtime = apt_runtime(), checks = length(checks),
  plan_hash = tools::md5sum(file.path(long_root, "inputs/plan.rds")), workers = workers), file.path(out, "audit.rds"))
writeLines("COMPLETE", file.path(out, "status.txt"))
cat("FRESH50K AUDIT PASSED; checks", length(checks), "; workers", workers,
  "; no real-data production chain or submission.\nOutput:", out, "\n")
