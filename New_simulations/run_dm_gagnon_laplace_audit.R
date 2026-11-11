args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

out_dir <- arg_value(
  "output",
  file.path("New_simulations", "outputs", "dm_gagnon_laplace_audit")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path("New_simulations", ".rcpp_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(
  file.path("New_simulations", "DM_RJMCMC_VS_sampler_online_adaptive.cpp"),
  cacheDir = cache_dir, rebuild = FALSE, showOutput = FALSE, verbose = FALSE
)

simulate_audit_data <- function(seed = 20291001L, n = 40L, J = 3L, P = 8L) {
  set.seed(seed)
  sigma <- outer(seq_len(P), seq_len(P), function(i, j) 0.55^abs(i - j))
  X <- scale(matrix(rnorm(n * P), n, P) %*% chol(sigma))
  beta <- matrix(0, P, J)
  beta[cbind(c(1, 4, 2, 6), c(1, 1, 2, 3))] <- c(0.6, -0.5, -0.6, 0.55)
  alpha <- log(c(1.3, 1.0, 0.8))
  gamma <- exp(sweep(X %*% beta, 2L, alpha, "+"))
  Y <- matrix(0L, n, J)
  for (i in seq_len(n)) {
    latent <- rgamma(J, shape = gamma[i, ], rate = 1)
    Y[i, ] <- as.integer(rmultinom(1L, 60L, latent / sum(latent)))
  }
  list(X = unclass(X), Y = Y, alpha = alpha, beta = beta)
}

dat <- simulate_audit_data()
P <- ncol(dat$X)
J <- ncol(dat$Y)
set.seed(20291002L)
states <- list(
  null = matrix(0L, P, J),
  sparse = matrix(rbinom(P * J, 1L, 0.15), P, J),
  moderate = matrix(rbinom(P * J, 1L, 0.35), P, J),
  dense = matrix(rbinom(P * J, 1L, 0.70), P, J)
)

rows <- vector("list", length(states))
for (state_index in seq_along(states)) {
  xi <- states[[state_index]]
  beta_full <- matrix(rnorm(P * J, sd = 0.20), P, J) * xi
  p25_a <- audit_dm_gagnon_laplace_probabilities(
    dat$X, dat$Y, dat$alpha, xi, beta_full, rep(10, J),
    a = 1, b = 9, optimizer_maxit = 25L
  )
  p25_b <- audit_dm_gagnon_laplace_probabilities(
    dat$X, dat$Y, dat$alpha, xi, beta_full, rep(10, J),
    a = 1, b = 9, optimizer_maxit = 25L
  )
  p100 <- audit_dm_gagnon_laplace_probabilities(
    dat$X, dat$Y, dat$alpha, xi, beta_full, rep(10, J),
    a = 1, b = 9, optimizer_maxit = 100L
  )
  rows[[state_index]] <- data.frame(
    state = names(states)[state_index],
    finite = all(is.finite(p25_a)) && all(is.finite(p100)),
    minimum_probability = min(p25_a),
    maximum_sum_error = max(abs(colSums(p25_a) - 1)),
    repeat_max_abs = max(abs(p25_a - p25_b)),
    reference_max_abs = max(abs(p25_a - p100))
  )
}

statewise <- do.call(rbind, rows)
tolerance_sum <- 1e-12
tolerance_repeat <- 1e-10
tolerance_reference <- 1e-2
summary <- data.frame(
  states = nrow(statewise),
  tolerance_sum = tolerance_sum,
  tolerance_repeat = tolerance_repeat,
  tolerance_reference = tolerance_reference,
  maximum_sum_error = max(statewise$maximum_sum_error),
  repeat_max_abs = max(statewise$repeat_max_abs),
  reference_max_abs = max(statewise$reference_max_abs),
  passed = all(statewise$finite) &&
    all(statewise$minimum_probability > 0) &&
    max(statewise$maximum_sum_error) <= tolerance_sum &&
    max(statewise$repeat_max_abs) <= tolerance_repeat &&
    max(statewise$reference_max_abs) <= tolerance_reference
)

write.csv(statewise, file.path(out_dir, "statewise_differences.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
writeLines(c(
  "Audit: Gagnon-adapted conditional Laplace/Barker model proposal",
  "Checks: finite positive probabilities, normalization, deterministic repetition,",
  "and agreement of 25-iteration production optimization with a 100-iteration reference.",
  "This is a numerical proposal audit; posterior recovery is assessed by the replicated experiment."
), file.path(out_dir, "manifest.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

print(summary)
if (!isTRUE(summary$passed[[1]])) {
  stop("Gagnon-Laplace proposal failed the numerical audit")
}
