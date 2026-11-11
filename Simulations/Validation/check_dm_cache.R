# Synthetic cache-equivalence smoke audit; run by the user, never on real data.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  x <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(x)) sub(paste0("^--", name, "="), "", x[1]) else default
}
out <- option("output", "Simulations/Validation/Results/New_Run/DM_Cache")
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a fresh output directory")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("Simulations/Validation/dm_validation.R")
Rcpp::sourceCpp("Simulations/Dirichlet_Multinomial/dm_sampler.cpp",
               cacheDir = "Simulations/.rcpp_cache", showOutput = FALSE)
set.seed(20261008L)
n <- 20L; P <- 4L; J <- 3L
X <- unclass(scale(matrix(rnorm(n * P), n, P)))
Y <- matrix(sample(0:12, n * J, replace = TRUE), n, J)
storage.mode(Y) <- "integer"
rows <- list()
near <- function(a, b) {
  if (!identical(dim(a), dim(b)) || length(a) != length(b)) return(FALSE)
  a <- as.numeric(a); b <- as.numeric(b)
  all(is.finite(a)) && all(is.finite(b)) &&
    max(abs(a - b), 0) <= 1e-10 * (1 + max(abs(b), 0))
}
for (budget in c(1L, 5L, 10L)) for (dense in c(FALSE, TRUE)) {
  initial <- matrix(as.integer(dense), P, J)
  run <- function(cache) {
    set.seed(20261009L + budget + 100L * dense)
    DM_online_adaptive_VS_RJMCMC(X, Y, rep(10, J), rep(10, J), 1, 9, 1, 1,
      runs = 120L, burn = 60L, thin = 1L, informed = TRUE, adapt_g = FALSE,
      initial_xi = initial, balance_directions = TRUE, uniform_mixture = .02,
      score_optimizer_maxit = budget, coefficient_optimizer_maxit = 3L,
      cache_score_fits = cache)
  }
  a <- run(FALSE); b <- run(TRUE)
  dm_validate_fit(a, "cache off"); dm_validate_fit(b, "cache on")
  for (field in c("xi", "alpha", "log_posterior_trace", "g_trace")) {
    rows[[length(rows) + 1L]] <- data.frame(budget = budget, dense = dense,
      check = field, passed = near(a[[field]], b[[field]]))
  }
  rows[[length(rows) + 1L]] <- data.frame(budget = budget, dense = dense,
    check = "active_coefficient_draws", passed = near(unlist(a$beta), unlist(b$beta)))
}
# Exercise exact-context fit reuse and cache invalidation without advancing a chain.
alpha <- c(.1, -.2, .3); context <- matrix(0, n, J)
design <- X[, 1:2, drop = FALSE]; start <- c(0, 0)
for (budget in c(1L, 5L, 10L)) {
  fit <- function(ctx) audit_dm_conditional_numerics(design, Y, alpha, ctx,
    1L, start, 2 / n, budget, TRUE, FALSE)
  first <- fit(context); repeat_fit <- fit(context)
  changed <- context; changed[, 2] <- X[, 3] * .2
  invisible(fit(changed)); after_history <- fit(context)
  for (field in c("mode", "mode_hessian")) {
    rows[[length(rows) + 1L]] <- data.frame(budget = budget, dense = NA,
      check = paste0("exact_context_", field),
      passed = near(first[[field]], repeat_fit[[field]]) && near(first[[field]], after_history[[field]]))
  }
}
checks <- do.call(rbind, rows)
write.csv(checks, file.path(out, "checks.csv"), row.names = FALSE)
write.csv(data.frame(passed = all(checks$passed)), file.path(out, "summary.csv"), row.names = FALSE)
hashes <- tools::md5sum(c(dm_source_files, "Simulations/Validation/check_dm_cache.R"))
write.csv(data.frame(file = names(hashes), md5 = unname(hashes)), file.path(out, "source_md5.csv"), row.names = FALSE)
print(checks)
if (!all(checks$passed)) stop("Cache equivalence audit failed")
cat("Cache audit passed. This is a finite synthetic equivalence check, not a general proof.\n")
