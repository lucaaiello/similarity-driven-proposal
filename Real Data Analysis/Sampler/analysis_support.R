source("Real Data Analysis/Sampler/runtime.R")
source("Real Data Analysis/Validation/Reference/Tempering/reference.R")
long_root <- "Real Data Analysis"
long_output <- "Real Data Analysis/Results/New_Run"
long_methods <- c("online_similarity_flip", "online_similarity_local", "uniform_flip")
long_names <- c("null", "moderate_5pct", "moderate_10pct", "dense_15pct")
long_files <- function() unique(c(aptopt_files(),apt_sources(),"Real Data Analysis/Sampler/analysis_support.R","Real Data Analysis/Sampler/ladder_steps.R","Real Data Analysis/prepare_inputs.R","Real Data Analysis/run_analysis.R","Real Data Analysis/summarize_results.R","Real Data Analysis/Validation/check_analysis.R"))
long_hashes <- function() tools::md5sum(long_files())
long_load <- function() {
  z<-readRDS("Real Data Analysis/Data/Prepared/analysis_plan.rds")
  apt_need(verify_recorded_sources(z$source_hashes) && verify_recorded_sources(z$original_hashes) &&
    identical(unname(z$input_hash),unname(tools::md5sum("Real Data Analysis/Data/Prepared/fresh_inputs.rds"))),
    "Recorded analysis plan/input/source receipt differs")
  list(plan=z$plan,inputs=readRDS("Real Data Analysis/Data/Prepared/fresh_inputs.rds"))
}
long_fresh <- function(x) {
  # Only X/Y, labels and the neighborhood are reused. No old physical state or tuning.
  x$origins <- list(); x$start_names <- long_names
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  P <- ncol(x$X); J <- ncol(x$Y)
  for (chain in 1:4) {
    set.seed(2026102900L + chain)
    xi <- matrix(0L, P, J)
    k <- as.integer(round(P * J * c(0, .05, .10, .15)[chain]))
    if (k) xi[sample.int(P * J, k)] <- 1L
    start <- list(alpha = rep(0, J), xi = xi, beta = lapply(colSums(xi), numeric),
      g = rep(1, J), Sigma_alpha = diag(J))
    for (method in long_methods) x$origins[[sprintf("%s_chain_%02d", method, chain)]] <-
      list(method = method, chain = chain, origin_start = long_names[chain], start = start,
        origin_iteration = 0L, recorded_log_posterior = NULL)
  }
  x
}
long_seed_scores <- function(x, method, initial) {
  # All ten replicas in a ladder begin at that ladder's genuinely fresh start.
  s <- x$origins[[sprintf("%s_chain_%02d", method, initial)]]$start
  ph <- apt_physical(s); eta <- matrix(0, nrow(x$Y), ncol(x$Y))
  ll <- get(".aptopt_cpp", .GlobalEnv)$audit_dm_stable_likelihood(x$Y, eta)
  lp <- sum(lgamma(ph$xi + 1) + lgamma(10 - ph$xi)) +
    sum(vapply(ph$beta, function(b) sum(dnorm(b, 0, sqrt(10), log = TRUE)), 0))
  v <- list(physical = ph, ll = ll, lp = lp, score = apt_score(x, ph))
  rep(list(v), 4L)
}
long_new <- function(x, plan, start, initial) {
  s <- apt_new_state(x, plan, start, FALSE, initial)
  s$schema <- "dm_apt50k_fresh_v1"; s$source_hashes <- long_hashes()
  s$input_hash <- tools::md5sum("Real Data Analysis/Data/Prepared/analysis_plan.rds")
  s$seed_origin[] <- start
  s$exchange_order <- matrix(NA_integer_, plan$iterations, plan$replicas - 1L)
  s$exchange_accepted <- matrix(FALSE, plan$iterations, plan$replicas - 1L)
  s
}
long_exchange <- function(s, x, iteration, beta) {
  L <- length(s$slots); phase <- if (iteration <= s$plan$burn) 1L else 2L
  ids <- match(c("3026", "4014"), rownames(x$X))
  # Random permutation of adjacent MH kernels. Recompute after EVERY attempt.
  # Every permutation and its reverse have equal selection probability.
  order <- sample.int(L - 1L)
  for (position in seq_along(order)) {
    edge <- order[position]
    before <- s$slots[[1]]$score$fitted[ids] > .05
    lr <- (beta[edge] - beta[edge + 1L]) *
      (s$slots[[edge + 1L]]$ll - s$slots[[edge]]$ll)
    accepted <- runif(1) < exp(min(0, lr))
    s$swaps_attempted[phase, edge] <- s$swaps_attempted[phase, edge] + 1L
    if (accepted) {
      s <- apt_swap_slots(s, edge)
      s$swaps_accepted[phase, edge] <- s$swaps_accepted[phase, edge] + 1L
      s$cold_swap_crossings[phase, ] <- s$cold_swap_crossings[phase, ] +
        as.integer(before != (s$slots[[1]]$score$fitted[ids] > .05))
    }
    s$exchange_order[iteration, position] <- edge
    s$exchange_accepted[iteration, position] <- accepted
    # Record endpoint arrivals within the exchange sweep, not only at its end.
    s$flow <- apt_flow_update(s$flow, s$labels, iteration)
  }
  s
}
long_gate <- function() {
  path <- "Real Data Analysis/Validation/Results/Analysis/audit.rds"
  apt_need(file.exists(path), "The matching cluster audit has not passed")
  a <- readRDS(path)
  apt_need(isTRUE(a$passed) && identical(a$source_hashes, long_hashes()) &&
    identical(a$runtime, apt_runtime()) && identical(a$plan_hash,
      tools::md5sum("Real Data Analysis/Data/Prepared/analysis_plan.rds")), "Fresh APT audit is failed/stale or from another runtime")
}
