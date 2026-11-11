args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

states <- as.integer(arg_value("states", "200"))
source("Simulations/Validation/dm_validation.R")
tolerance <- as.numeric(arg_value("tolerance", "1e-12"))
out_dir <- arg_value(
  "output", file.path("Simulations/Validation/Results/New_Run/DM_Adaptation")
)
if (states < 3L || !is.finite(tolerance) || tolerance <= 0) stop("invalid audit settings")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path("Simulations/.rcpp_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(
  file.path("Simulations/Dirichlet_Multinomial/dm_sampler.cpp"),
  cacheDir = cache_dir, rebuild = FALSE, showOutput = FALSE, verbose = FALSE
)

reference_probabilities <- function(log_weights, indicators, balanced, mixture) {
  normalize_group <- function(index, mass) {
    if (!length(index)) return(numeric())
    z <- exp(log_weights[index] - max(log_weights[index]))
    mass * z / sum(z)
  }
  p <- numeric(length(indicators))
  births <- which(indicators == 0L)
  deaths <- which(indicators == 1L)
  if (!balanced) {
    z <- exp(log_weights - max(log_weights))
    p <- z / sum(z)
  } else {
    birth_mass <- if (length(births) && length(deaths)) 0.5 else as.numeric(length(births) > 0L)
    death_mass <- if (length(births) && length(deaths)) 0.5 else as.numeric(length(deaths) > 0L)
    p[births] <- normalize_group(births, birth_mass)
    p[deaths] <- normalize_group(deaths, death_mass)
  }
  if (mixture > 0) {
    p <- (1 - mixture) * p
    if (balanced) {
      birth_mass <- if (length(births) && length(deaths)) 0.5 else as.numeric(length(births) > 0L)
      death_mass <- if (length(births) && length(deaths)) 0.5 else as.numeric(length(deaths) > 0L)
      if (length(births)) p[births] <- p[births] + birth_mass * mixture / length(births)
      if (length(deaths)) p[deaths] <- p[deaths] + death_mass * mixture / length(deaths)
    } else {
      p <- p + mixture / length(p)
    }
  }
  p
}

set.seed(20282001L)
max_abs <- 0
max_sum_error <- 0
max_direction_error <- 0
legacy_max_abs <- 0
for (state in seq_len(states)) {
  P <- sample(3:50, 1L)
  indicators <- if (state == 1L) rep(0L, P) else if (state == 2L) rep(1L, P) else rbinom(P, 1L, 0.25)
  log_weights <- rnorm(P, sd = sample(c(0.1, 1, 10), 1L))
  mixture <- sample(c(0, 0.02, 0.10), 1L)
  cpp <- audit_dm_similarity_normalization(log_weights, indicators, TRUE, mixture)
  reference <- reference_probabilities(log_weights, indicators, TRUE, mixture)
  max_abs <- max(max_abs, max(abs(cpp - reference)))
  max_sum_error <- max(max_sum_error, abs(sum(cpp) - 1))
  if (any(indicators == 0L) && any(indicators == 1L)) {
    max_direction_error <- max(
      max_direction_error,
      abs(sum(cpp[indicators == 0L]) - 0.5),
      abs(sum(cpp[indicators == 1L]) - 0.5)
    )
  }
  cpp_legacy <- audit_dm_similarity_normalization(log_weights, indicators, FALSE, 0)
  reference_legacy <- reference_probabilities(log_weights, indicators, FALSE, 0)
  legacy_max_abs <- max(legacy_max_abs, max(abs(cpp_legacy - reference_legacy)))
}

simulate_smoke_data <- function() {
  set.seed(20282002L)
  n <- 20L
  P <- 4L
  J <- 3L
  X <- scale(matrix(rnorm(n * P), n, P))
  beta <- matrix(0, P, J)
  beta[1, 1] <- 0.5
  beta[2, 2] <- -0.5
  eta <- X %*% beta
  gamma <- exp(eta)
  Y <- matrix(0L, n, J)
  for (i in seq_len(n)) {
    latent <- rgamma(J, gamma[i, ], 1)
    Y[i, ] <- as.integer(rmultinom(1L, 30L, latent / sum(latent)))
  }
  list(X = unclass(X), Y = Y)
}

dat <- simulate_smoke_data()
set.seed(20282003L)
smoke <- DM_online_adaptive_VS_RJMCMC(
  dat$X, dat$Y, s2 = rep(10, 3), r2 = rep(10, 3),
  a = 1, b = 9, g_initial = 1, c = 1,
  runs = 60L, burn = 30L, thin = 1L,
  informed = TRUE, adapt_g = TRUE, adaptation_start = 1L,
  adaptation_end = -1L, window_size = 5L, rm_scale = 0.5,
  decay_power = 0.6, min_g = 0.05, max_g = 10,
  use_move = FALSE, adaptation_metric = "diversity",
  target_effective_neighbors = 3, adaptation_offset = 10,
  balance_directions = TRUE, uniform_mixture = 0.02,
  score_optimizer_maxit = 1L, coefficient_optimizer_maxit = 2L,
  cache_score_fits = TRUE
)

smoke_checks <- c(
  finite_posterior = all(is.finite(smoke$log_posterior_trace)),
  finite_alpha = all(is.finite(smoke$alpha)),
  finite_beta = all(is.finite(unlist(smoke$beta))),
  finite_g = all(is.finite(smoke$g_trace)),
  bounded_g = all(smoke$g_trace >= 0.05 & smoke$g_trace <= 10),
  finite_effective_neighbors = all(is.finite(smoke$proposal_effective_neighbors_trace)),
  valid_effective_neighbors = all(smoke$proposal_effective_neighbors_trace >= 1 &
                                    smoke$proposal_effective_neighbors_trace <= ncol(dat$X) + tolerance),
  adaptation_updated = all(smoke$adaptation_updates_by_category > 0),
  retained_dimensions = identical(dim(smoke$xi), c(ncol(dat$X), ncol(dat$Y), 30L))
)

passed <- max_abs <= tolerance && legacy_max_abs <= tolerance &&
  max_sum_error <= tolerance && max_direction_error <= tolerance && all(smoke_checks)
summary <- data.frame(
  states = states, tolerance = tolerance, max_abs = max_abs,
  legacy_max_abs = legacy_max_abs, max_sum_error = max_sum_error,
  max_direction_error = max_direction_error,
  smoke_checks_passed = all(smoke_checks), passed = passed
)
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(data.frame(check = names(smoke_checks), passed = unname(smoke_checks)),
          file.path(out_dir, "smoke_checks.csv"), row.names = FALSE)
writeLines(c(
  "Audit scope:",
  "- independent normalization check for balanced and legacy one-flip probabilities",
  "- direction-mass and probability-sum invariants",
  "- short online-adaptation smoke test and output-shape checks",
  "- this is not a Gaussian rank-one audit; the DM PMLE/Hessian calculations remain numerical"
), file.path(out_dir, "manifest.txt"))
print(summary)
if (!passed) stop("DM online-adaptation audit failed")
