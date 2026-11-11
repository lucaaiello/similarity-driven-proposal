source("Real Data Analysis/Provenance/verify_sources.R")
apt_root <- "Real Data Analysis/Validation/Reference/Tempering"
apt_need <- function(ok, message) { if (!isTRUE(ok)) stop(message, call. = FALSE); invisible(TRUE) }
apt_options <- function(allowed) {
  a <- commandArgs(TRUE); k <- sub("^--([^=]+)=.*$", "\\1", a)
  apt_need(all(grepl("^--[^=]+=.+$", a)) && !anyDuplicated(k) && all(k %in% allowed), "Use unique --name=value arguments")
  setNames(sub("^--[^=]+=", "", a), k)
}
apt_get <- function(o, k, default) if (k %in% names(o)) o[[k]] else default
apt_environment <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1", OMP_DYNAMIC = "FALSE",
    MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
}
apt_runtime <- function() list(R = R.version.string, platform = R.version$platform,
  packages = vapply(c("Rcpp", "RcppArmadillo", "RcppEnsmallen"), function(p) as.character(packageVersion(p)), ""),
  threads = as.integer(Sys.getenv("OMP_NUM_THREADS", "1")))
apt_sources <- function() file.path(apt_root,c("dm_sampler_reference.cpp","joint_kernel.h","dm_stable_numerics.h","support.R","math.R","runner.R","reference.R","Untempered/dm_sampler_untempered.cpp","Untempered/joint_kernel.h","Untempered/dm_stable_numerics.h"))
apt_hashes <- function() { p <- apt_sources(); apt_need(all(file.exists(p)), "Missing APT source"); tools::md5sum(p) }
apt_load <- function() {
  z<-readRDS("Real Data Analysis/Data/Reference/reference_plan.rds")
  apt_need(verify_recorded_sources(z$source_hashes) && identical(unname(z$origin_hash),
    unname(tools::md5sum("Real Data Analysis/Data/Reference/reference_states.rds"))),"Original reference input/source receipt differs")
  list(plan=z$plan,inputs=readRDS("Real Data Analysis/Data/Reference/reference_states.rds"),receipt=z)
}
apt_compile <- function(baseline = FALSE) {
  e <- new.env(parent = globalenv())
  Rcpp::sourceCpp(file.path(apt_root, if (baseline) "Untempered/dm_sampler_untempered.cpp" else "dm_sampler_reference.cpp"),
    env = e, cacheDir = file.path(tempdir(), paste0("apt_", Sys.getpid()), if (baseline) "baseline" else "new"), showOutput = FALSE)
  e$dm_runtime_use_cached(TRUE); e$dm_runtime_v2_options(TRUE, TRUE, TRUE); e$dm_runtime_profile_reset(FALSE)
  e
}
apt_args <- function(x, physical, tuning, inverse_temperature, sweeps = 1L, method = "online_similarity_local") {
  informed <- method != "uniform_flip"; local <- method == "online_similarity_local"
  # Fixed temperature-dependent proposal rule, leaving all priors unchanged.
  start <- c(physical, list(g = tuning$g,
    Sigma_alpha = tuning$Sigma_alpha * min(1 / inverse_temperature, 100)))
  list(X = x$X, Y = x$Y, s2 = rep(10, ncol(x$Y)), r2 = rep(10, ncol(x$Y)),
    a = 1, b = 9, g_initial = 1, c = 1, runs = as.integer(sweeps), burn = 0L, thin = 1L,
    informed = informed, adapt_g = FALSE, adaptation_start = 1L, adaptation_end = 0L,
    balance_directions = informed, uniform_mixture = if (informed) .02 else 0,
    neigh = if (local) x$neighbors else NULL, use_move = local, n_move_attempts = 1L,
    score_optimizer_maxit = 1L, coefficient_optimizer_maxit = 25L,
    initial_xi = physical$xi, cache_score_fits = TRUE, coefficient_warm_start = TRUE,
    within_model_update = TRUE, within_model_attempts = 1L, within_model_scale = .25,
    within_model_fixed_probability = .5, within_model_fixed_scale = .05,
    joint_update = TRUE, joint_scale = 2, joint_attempts = 1L, pilot_start = start,
    inverse_temperature = inverse_temperature, quiet = TRUE)
}
apt_physical <- function(s) s[c("alpha", "xi", "beta")]
apt_score <- function(x, physical) {
  P <- ncol(x$X); J <- ncol(x$Y); B <- matrix(0, P, J)
  for (j in seq_len(J)) B[physical$xi[, j] == 1L, j] <- physical$beta[[j]]
  eta <- sweep(x$X %*% B, 2, physical$alpha, "+")
  top <- apply(eta, 1, max); total <- top + log(rowSums(exp(sweep(eta, 1, top, "-"))))
  jp <- match("Prevotella", sub("^.*\\.", "", x$taxa))
  fitted <- exp(eta[, jp] - total)
  list(prev_alpha = physical$alpha[jp], prev_size = sum(physical$xi[, jp]),
       model_size = sum(physical$xi), fitted = fitted, mean_log_total = mean(total))
}
apt_worker_initialize <- function(x) {
  apt_environment()
  assign(".apt_x", x, envir = .GlobalEnv)
  assign(".apt_cpp", apt_compile(FALSE), envir = .GlobalEnv)
  list(runtime = apt_runtime(), pid = Sys.getpid())
}
apt_worker_initial_scores <- function(method) {
  x <- get(".apt_x", envir = .GlobalEnv); e <- get(".apt_cpp", envir = .GlobalEnv)
  lapply(1:4, function(chain) {
    z <- x$origins[[sprintf("%s_chain_%02d", method, chain)]]
    physical <- apt_physical(z$start)
    B <- matrix(0, ncol(x$X), ncol(x$Y))
    for (j in seq_len(ncol(x$Y))) B[physical$xi[, j] == 1L, j] <- physical$beta[[j]]
    ll <- e$audit_dm_stable_likelihood(x$Y, sweep(x$X %*% B, 2, physical$alpha, "+"))
    lp <- -.5 * sum(physical$alpha^2 / 10) + sum(vapply(physical$beta,
      function(b) sum(dnorm(b, 0, sqrt(10), log = TRUE)), 0)) +
      sum(lgamma(physical$xi + 1) + lgamma(10 - physical$xi))
    if (!is.null(z$recorded_log_posterior)) apt_need(abs(ll + lp - z$recorded_log_posterior) < 1e-4,
      "Initial tempered state differs from saved posterior target")
    list(physical = physical, ll = ll, lp = lp, score = apt_score(x, physical))
  })
}
apt_worker_step <- function(request) {
  cpu_before <- sum(proc.time()[1:2])
  assign(".Random.seed", request$rng, envir = .GlobalEnv)
  x <- get(".apt_x", envir = .GlobalEnv); e <- get(".apt_cpp", envir = .GlobalEnv)
  args <- apt_args(x, request$physical, request$tuning, request$beta, method = request$method)
  f <- do.call(e$DM_tempered_VS_RJMCMC, args)
  cs <- f$continuation_state; physical <- apt_physical(cs)
  apt_need(f$nonfinite_mh_ratios_rejected == 0 && f$optimizer_zero_fallbacks == 0 &&
    is.finite(cs$log_lik) && all(is.finite(physical$alpha)) &&
    all(is.finite(unlist(physical$beta))) && all(physical$xi %in% 0:1), "Nonfinite/invalid local APT update")
  list(physical = physical, ll = cs$log_lik,
    lp = cs$log_prior_alpha + sum(cs$log_prior_beta) + sum(cs$log_prior_xi),
    rng = get(".Random.seed", envir = .GlobalEnv),
    effective = as.numeric(f$proposal_effective_neighbors_trace[1, ]),
    counts = c(alpha = f$acpt_rate_alpha, flip = sum(f$acpt_count_betaxi_by_category),
      local = sum(f$acpt_count_move_by_category), joint = sum(f$joint_accepted_by_category),
      projections = f$precision_projection_calls),
    local_seconds = f$total_seconds, score = apt_score(x, physical),
    cpu_seconds = sum(proc.time()[1:2]) - cpu_before)
}
apt_pool <- function(x, workers) {
  apt_need(workers >= 1L && workers <= 16L, "Invalid worker count")
  cl <- parallel::makePSOCKcluster(workers, outfile = "")
  ok <- FALSE
  on.exit(if (!ok) parallel::stopCluster(cl))
  parallel::clusterCall(cl, function(wd) {
    setwd(wd); source("Real Data Analysis/Validation/Reference/Tempering/support.R")
    source("Real Data Analysis/Validation/Reference/Tempering/math.R"); NULL
  }, getwd())
  runtime <- parallel::clusterCall(cl, apt_worker_initialize, x)
  apt_need(all(vapply(runtime, function(v) identical(v$runtime, apt_runtime()), TRUE)), "Worker runtimes differ")
  ok <- TRUE; attr(cl, "runtime") <- runtime; cl
}
apt_write <- function(z, path) {
  apt_need(!file.exists(path), paste("Refusing overwrite:", path))
  temp <- paste0(path, ".tmp_", Sys.getpid())
  apt_need(!file.exists(temp) && dir.create(temp, recursive = TRUE), "Cannot create fresh checkpoint staging")
  saveRDS(z, file.path(temp, "state.rds"))
  writeLines(unname(tools::md5sum(file.path(temp, "state.rds"))), file.path(temp, "state.md5"))
  writeLines("COMPLETE", file.path(temp, "status.txt"))
  apt_need(file.rename(temp, path), "Checkpoint install failed; preserve staging")
}
apt_read <- function(path) {
  p <- file.path(path, "state.rds")
  apt_need(identical(readLines(file.path(path, "status.txt")), "COMPLETE") &&
    identical(readLines(file.path(path, "state.md5")), unname(tools::md5sum(p))), "Invalid APT checkpoint")
  readRDS(p)
}
apt_gate <- function() {
  p <- "Real Data Analysis/Validation/Results/Reference/audit.rds"; apt_need(file.exists(p), "Full cluster APT audit is missing")
  a <- readRDS(p)
  apt_need(isTRUE(a$passed) && identical(a$source_hashes, apt_hashes()) &&
    identical(a$input_hash, tools::md5sum("Real Data Analysis/Data/Reference/reference_plan.rds")) &&
    identical(a$runtime, apt_runtime()), "APT audit failed, stale or from a different runtime")
}
