args <- commandArgs(trailingOnly = TRUE)
hit <- grep("^--output=", args, value = TRUE)
out_dir <- if (length(hit)) sub("^--output=", "", hit[[1]]) else
  "Simulations/Validation/Results/New_Run/DM_Numerics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("numerical audit: loading validation helpers")
source("Simulations/Validation/dm_validation.R")
# A failed rerun must never leave an earlier passed gate in place.
write.csv(data.frame(passed = FALSE, numerical_engine = dm_numerical_engine),
          file.path(out_dir, "summary.csv"), row.names = FALSE)
message("numerical audit: compiling/loading C++ module")
Rcpp::sourceCpp("Simulations/Dirichlet_Multinomial/dm_sampler.cpp",
                cacheDir = "Simulations/.rcpp_cache", rebuild = FALSE)
message("numerical audit: C++ module loaded")
checks <- list()
check <- function(name, error = 0, tolerance = 0, valid = TRUE) {
  passed <- isTRUE(valid) && is.finite(error) && error <= tolerance
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, error = error, tolerance = tolerance, passed = passed
  )
  if (!passed) stop("Numerical audit failed: ", name, " (error=", error, ")")
}
relative_error <- function(x, y) max(abs(x - y) / (1 + abs(y)))
logadd <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))
# Independent integer-count reference, without gamma/digamma subtraction.
rising_reference <- function(log_a, m) {
  if (!m) return(c(0, 0, 0))
  k <- seq_len(m) - 1
  log_k <- log(k)
  ratio <- plogis(log_a - log_k)
  c(sum(logadd(log_a, log_k) - log_a),
    sum(ratio), sum(ratio * (1 - ratio)))
}
cases <- expand.grid(log_a = c(-1000, -40, -36, -10, 0, 10, 15, 20, 40, 1000),
                     count = c(0L, 1L, 2L, 64L, 65L, 1000L, 65000L))
cases <- rbind(cases, data.frame(
  log_a = log(65000 * 100) + c(-1e-6, 1e-6), count = 65000L
))
reference <- t(vapply(seq_len(nrow(cases)), function(i)
  rising_reference(cases$log_a[i], cases$count[i]), numeric(3)))
message("numerical audit: rising-factorial reference")
cpp <- audit_dm_rising_values(cases$log_a, as.integer(cases$count))
check("rising_factorial_and_derivatives", relative_error(cpp, reference), 2e-7)

likelihood_reference <- function(Y, eta) {
  sum(vapply(seq_len(nrow(Y)), function(i) {
    ls <- max(eta[i, ]) + log(sum(exp(eta[i, ] - max(eta[i, ]))))
    own <- vapply(seq_len(ncol(Y)), function(j) {
      m <- Y[i, j]
      if (!m) return(0)
      sum(logadd(eta[i, j], log(seq_len(m) - 1)))
    }, numeric(1))
    total <- sum(Y[i, ])
    sum(own) - if (total) sum(logadd(ls, log(seq_len(total) - 1))) else 0
  }, numeric(1)))
}
set.seed(20261005L)
message("numerical audit: likelihood references")
n <- 12L; P <- 5L; J <- 3L
X <- unclass(scale(matrix(rnorm(n * P), n, P)))
Y <- matrix(sample(0:100, n * J, replace = TRUE), n, J)
storage.mode(Y) <- "integer"
alpha <- c(0.2, -0.3, 0.4)
xbeta <- matrix(rnorm(n * J, sd = 0.1), n, J)
for (shift in c(-1000, -40, 0, 40, 1000)) {
  eta <- sweep(xbeta, 2L, alpha + shift, "+")
  ref <- likelihood_reference(Y, eta)
  check(paste0("likelihood_shift_", shift),
        relative_error(audit_dm_stable_likelihood(Y, eta), ref), 2e-8)
}
eta <- sweep(xbeta, 2L, alpha, "+")
gamma <- exp(eta); sums <- rowSums(gamma); counts <- rowSums(Y)
legacy <- sum(lgamma(sums) - lgamma(sums + counts)) +
  sum(lgamma(Y + gamma) - lgamma(gamma))
check("moderate_legacy_likelihood", abs(audit_dm_stable_likelihood(Y, eta) - legacy), 1e-8)

beta <- c(0.2, -0.1); design <- X[, 1:2]; penalty <- 0.2
message("numerical audit: conditional derivatives")
conditional <- function(b, offsets = alpha) audit_dm_conditional_numerics(
  design, Y, offsets, xbeta, 1L, b, penalty, 3L, FALSE
)
for (shift in c(-40, 0, 40, 1000)) {
  message("numerical audit: conditional shift ", shift, " (objective)")
  offsets <- alpha + shift
  obj <- conditional(beta, offsets)
  ref_eta <- sweep(xbeta, 2L, offsets, "+")
  old <- ref_eta
  ref_eta[, 1] <- offsets[1] + drop(design %*% beta)
  ref_value <- -likelihood_reference(Y, ref_eta) + penalty * sum(beta^2) / 2
  check(paste0("objective_", shift), relative_error(obj$value, ref_value), 2e-8)
  check(paste0("incremental_delta_", shift), relative_error(obj$delta,
    likelihood_reference(Y, ref_eta) - likelihood_reference(Y, old)), 2e-7)
  message("numerical audit: conditional shift ", shift, " (finite differences)")
  h <- 1e-4
  numeric_gradient <- numeric(2); numeric_hessian <- matrix(0, 2, 2)
  for (k in 1:2) {
    plus <- minus <- beta; plus[k] <- plus[k] + h; minus[k] <- minus[k] - h
    p <- conditional(plus, offsets); m <- conditional(minus, offsets)
    numeric_gradient[k] <- (p$value - m$value) / (2 * h)
    numeric_hessian[, k] <- -(p$gradient - m$gradient) / (2 * h)
  }
  check(paste0("gradient_", shift), relative_error(obj$gradient, numeric_gradient), 5e-5)
  check(paste0("hessian_", shift), relative_error(obj$hessian, numeric_hessian), 5e-5)
}
# The optimizer is intentionally tested only in the ordinary finite regime.
# Extreme synthetic concentrations above test the target and its derivatives,
# not a finite-iteration quasi-Newton line search.
optimized_a <- audit_dm_conditional_numerics(design, Y, alpha, xbeta, 1L,
  beta, penalty, 3L, TRUE)
optimized_b <- audit_dm_conditional_numerics(design, Y, alpha, xbeta, 1L,
  beta + 20, penalty, 3L, TRUE)
check("moderate_optimizer_history_independence",
      max(abs(optimized_a$mode - optimized_b$mode)), 1e-12)

message("numerical audit: Gaussian proposal consistency")
for (type in c("spd", "indefinite", "asymmetric")) {
  Q <- switch(type, spd = matrix(c(2, 0.3, 0.3, 1), 2),
              indefinite = matrix(c(1, 2, 2, 1), 2),
              asymmetric = matrix(c(2, 0.31, 0.3, 1), 2))
  mean <- c(0.2, -0.1); x <- c(0.5, 0.8)
  result <- audit_dm_gaussian_numerics(mean, Q, x, 20000L)
  stabilized <- result$precision
  ref <- -log(2 * pi) + as.numeric(determinant(stabilized, logarithm = TRUE)$modulus) / 2 -
    drop(t(x - mean) %*% stabilized %*% (x - mean)) / 2
  check(paste0("gaussian_density_", type), abs(result$log_density - ref), 1e-10)
  whitening <- chol(stabilized) %*% sweep(result$samples, 1L, mean, "-")
  check(paste0("gaussian_sample_covariance_", type), max(abs(cov(t(whitening)) - diag(2))), 0.06)
  check(paste0("gaussian_sample_mean_", type), max(abs(rowMeans(whitening))), 0.06)
}
bad_precision <- try(audit_dm_gaussian_numerics(c(0, 0), diag(c(Inf, 1)), c(0, 0)), silent = TRUE)
check("nonfinite_precision_fails_closed", valid = inherits(bad_precision, "try-error"))
accept <- audit_dm_safe_acceptance(c(NaN, Inf, -Inf, 0, -2, 3))
check("nonfinite_mh_rejected", valid = identical(as.numeric(accept), c(-Inf, -Inf, -Inf, 0, -2, 0)))

adjacency <- matrix(1L, P, P); diag(adjacency) <- 0L
message("numerical audit: flip/local proposal references")
for (pattern in c("null", "mixed", "full")) {
  xi <- matrix(0L, P, J)
  xi[, 1] <- switch(pattern, null = rep(0L, P), mixed = c(1L, 1L, 0L, 0L, 0L), full = rep(1L, P))
  b <- matrix(0, P, J); b[, 1] <- rnorm(P, sd = 0.1) * xi[, 1]
  probs <- audit_dm_neighborhood_probabilities(X, Y, alpha, xi, b, 1L, adjacency)
  # Independent weights, using production fitted points but an independent
  # integer-count likelihood and R's pchisq/softmax (not C++ normalization).
  score <- function(indicators) {
    selected <- which(indicators == 1L)
    linear <- rep(0, n)
    if (length(selected)) {
      fit <- audit_dm_conditional_numerics(X[, selected, drop = FALSE], Y, alpha,
        matrix(0, n, J), 1L, rep(0, length(selected)), length(selected) / n, 1L, TRUE)
      linear <- drop(X[, selected, drop = FALSE] %*% fit$mode)
    }
    e <- matrix(rep(alpha, each = n), n, J); base <- likelihood_reference(Y, e)
    e[, 1] <- e[, 1] + linear
    stat <- max(0, 2 * (likelihood_reference(Y, e) - base))
    if (!length(selected)) return(0)
    -pchisq(stat, length(selected), lower.tail = FALSE, log.p = TRUE) / log(10)
  }
  softmax <- function(z) { e <- exp(z - max(z)); e / sum(e) }
  scores <- vapply(seq_len(P), function(p) {
    proposal <- xi[, 1]; proposal[p] <- 1L - proposal[p]; score(proposal)
  }, numeric(1))
  ref_flip <- numeric(P)
  groups <- Filter(length, list(which(xi[, 1] == 0L), which(xi[, 1] == 1L)))
  for (group in groups) ref_flip[group] <-
    (0.98 * softmax(scores[group]) + 0.02 / length(group)) / length(groups)
  check(paste0("flip_probability_reference_", pattern), max(abs(probs$flip - ref_flip)), 2e-7)
  ref_swap <- matrix(0, P, P)
  active <- which(xi[, 1] == 1L); inactive <- which(xi[, 1] == 0L)
  if (length(active) && length(inactive)) for (drop in active) {
    z <- vapply(inactive, function(add) {
      proposal <- xi[, 1]; proposal[drop] <- 0L; proposal[add] <- 1L; score(proposal)
    }, numeric(1))
    ref_swap[drop, inactive] <- softmax(z) / length(active)
  }
  check(paste0("local_probability_reference_", pattern), max(abs(probs$swap - ref_swap)), 2e-7)
}

# Exercise all three actual kernels, including local swaps, on count data.
message("numerical audit: finite-state sampler smoke checks")
for (method in c("flip", "uniform", "local")) for (start in c("null", "dense")) {
  set.seed(20261006L + match(method, c("flip", "uniform", "local")) * 10L + (start == "dense"))
  initial <- if (start == "null") matrix(0L, P, J) else matrix(rbinom(P * J, 1, 0.5), P, J)
  fit <- DM_online_adaptive_VS_RJMCMC(X, Y, rep(10, J), rep(10, J), 1, 9, 1, 1,
    runs = 120L, burn = 60L, informed = method != "uniform", adapt_g = method != "uniform",
    adaptation_start = 1L, adaptation_end = 60L, window_size = 10L,
    adaptation_metric = "diversity", target_effective_neighbors = 3,
    balance_directions = method != "uniform", uniform_mixture = if (method == "uniform") 0 else 0.02,
    neigh = lapply(seq_len(P), function(p) setdiff(seq_len(P), p)), use_move = method == "local",
    initial_xi = initial, score_optimizer_maxit = 1L, coefficient_optimizer_maxit = 3L)
  dm_validate_fit(fit, paste(method, start))
  check(paste0("finite_sampler_", method, "_", start))
  check(paste0("no_nonfinite_mh_ratios_", method, "_", start),
        valid = fit$nonfinite_mh_ratios_rejected == 0L)
  if (method == "local") check(paste0("swap_exercised_", start), valid = sum(fit$prop_count_move_by_category) > 0)
}
synthetic <- fit; synthetic$log_posterior_trace[1] <- Inf
check("invalid_checkpoint_rejected", valid = inherits(try(dm_validate_fit(synthetic), silent = TRUE), "try-error"))
table <- do.call(rbind, checks)
write.csv(table, file.path(out_dir, "checks.csv"), row.names = FALSE)
hashes <- tools::md5sum(dm_source_files)
write.csv(data.frame(file = names(hashes), md5 = unname(hashes)),
          file.path(out_dir, "source_md5.csv"), row.names = FALSE)
summary <- data.frame(checks = nrow(table), passed = all(table$passed),
  numerical_engine = dm_numerical_engine, fit_engine = dm_fit_engine)
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
writeLines(c("Independent integer-count likelihood/derivative references, finite differences,",
  "Gaussian sampling/density consistency, statewise flip/local probability references,",
  "and all-kernel finite-state smoke checks. Not a real-data convergence proof or a",
  "full trans-dimensional detailed-balance proof."), file.path(out_dir, "scope.txt"))
print(summary)
