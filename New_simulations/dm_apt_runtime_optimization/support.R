# Verification/integration only. Existing sources and scientific results are immutable.
source("New_simulations/dm_apt_runtime_optimization/baseline/support.R")
source("New_simulations/dm_apt_runtime_optimization/baseline/math.R")
source("New_simulations/dm_apt_runtime_optimization/baseline/runner.R")
source("New_simulations/dm_apt_runtime_optimization/audit_cases.R")
aptopt_root <- "New_simulations/dm_apt_runtime_optimization"
aptopt_protected <- function() {
  # Original complete APT source graph, inputs and completed checkpoint receipts.
  p <- c(apt_sources(), file.path(apt_root, "inputs/plan.rds"), file.path(apt_root, "inputs/origins.rds"))
  summary <- "New_simulations/outputs/dm_apt10_v1/summary/summary.rds"
  if (file.exists(summary)) {
    s <- readRDS(summary); p <- c(p, summary, names(s$checkpoint_hashes))
  }
  tools::md5sum(p)
}
aptopt_files <- function() file.path(aptopt_root, c("DM_persistent.cpp", "joint_kernel.h", "dm_stable_numerics.h",
  "baseline/DM_tempered.cpp", "baseline/joint_kernel.h", "baseline/dm_stable_numerics.h",
  "baseline/support.R", "baseline/math.R", "baseline/runner.R", "baseline/reference.R",
  "support.R", "audit_cases.R", "verify.R"))
aptopt_hashes <- function() tools::md5sum(aptopt_files())
aptopt_compile <- function(engine) {
  apt_need(engine %in% c("reference", "persistent"), "Unknown optimization engine")
  path <- if (engine == "reference") "baseline/DM_tempered.cpp" else "DM_persistent.cpp"
  e <- new.env(parent = globalenv())
  Rcpp::sourceCpp(file.path(aptopt_root, path), env = e, showOutput = FALSE,
    cacheDir = file.path(tempdir(), paste0("aptopt_", Sys.getpid()), engine))
  e$dm_runtime_use_cached(TRUE); e$dm_runtime_v2_options(TRUE, TRUE, TRUE); e$dm_runtime_profile_reset(FALSE)
  e
}
aptopt_initialize <- function(x, engine = "both") {
  apt_environment()
  assign(".apt_x", x, envir = .GlobalEnv)
  if (engine %in% c("both", "reference")) assign(".apt_cpp", aptopt_compile("reference"), .GlobalEnv)
  if (engine %in% c("both", "persistent")) {
    assign(".aptopt_cpp", aptopt_compile("persistent"), .GlobalEnv)
    assign(".aptopt_contexts", new.env(parent = emptyenv()), .GlobalEnv)
  }
  assign(".aptopt_engine", engine, .GlobalEnv)
  list(runtime = apt_runtime(), pid = Sys.getpid())
}
aptopt_context <- function(method) {
  contexts <- get(".aptopt_contexts", .GlobalEnv)
  if (!exists(method, contexts, inherits = FALSE)) {
    x <- get(".apt_x", .GlobalEnv); e <- get(".aptopt_cpp", .GlobalEnv)
    assign(method, e$apt_context_create(x$X, x$Y, x$neighbors, method), contexts)
  }
  get(method, contexts, inherits = FALSE)
}
aptopt_step <- function(request) {
  engine <- request$engine
  if (is.null(engine)) engine <- get(".aptopt_engine", .GlobalEnv)
  apt_need(engine %in% c("reference", "persistent"), "Specify reference or persistent engine")
  cap <- if (is.null(request$coefficient_cap)) 25L else request$coefficient_cap
  score_cap <- if (is.null(request$score_cap)) 1L else request$score_cap
  if (engine == "reference" && cap == 25L && score_cap == 1L) return(apt_worker_step(request))
  cpu_before <- sum(proc.time()[1:2]); assign(".Random.seed", request$rng, .GlobalEnv)
  x <- get(".apt_x", .GlobalEnv)
  if (engine == "reference") {
    args <- apt_args(x, request$physical, request$tuning, request$beta, method = request$method)
    args$coefficient_optimizer_maxit <- cap; args$score_optimizer_maxit <- score_cap
    f <- do.call(get(".apt_cpp", .GlobalEnv)$DM_tempered_VS_RJMCMC, args)
    cs <- f$continuation_state; ph <- apt_physical(cs)
    apt_need(f$nonfinite_mh_ratios_rejected == 0 && f$optimizer_zero_fallbacks == 0, "Reference numerical failure")
    answer <- list(physical = ph, ll = cs$log_lik,
      lp = cs$log_prior_alpha + sum(cs$log_prior_beta) + sum(cs$log_prior_xi),
      rng = get(".Random.seed", .GlobalEnv), effective = as.numeric(f$proposal_effective_neighbors_trace[1, ]),
      counts = c(alpha = f$acpt_rate_alpha, flip = sum(f$acpt_count_betaxi_by_category),
        local = sum(f$acpt_count_move_by_category), joint = sum(f$joint_accepted_by_category), projections = f$precision_projection_calls),
      local_seconds = f$total_seconds, score = apt_score(x, ph))
  } else {
    e <- get(".aptopt_cpp", .GlobalEnv)
    f <- e$apt_context_step(aptopt_context(request$method), request$physical, request$tuning$g,
      request$tuning$Sigma_alpha * min(1 / request$beta, 100), request$beta,
      as.integer(cap), as.integer(score_cap), 2)
    apt_need(f$nonfinite_mh_ratios_rejected == 0 && f$optimizer_zero_fallbacks == 0 &&
      all(is.finite(f$physical$alpha)) && all(is.finite(unlist(f$physical$beta))) &&
      all(f$physical$xi %in% 0:1) && is.finite(f$ll) && is.finite(f$lp), "Persistent numerical failure")
    # Match the original R summation order exactly, independently of native lp.
    lp <- f$prior_alpha + sum(f$prior_beta) + sum(f$prior_xi)
    answer <- list(physical = f$physical, ll = f$ll, lp = lp,
      rng = get(".Random.seed", .GlobalEnv), effective = as.numeric(f$effective),
      counts = f$counts, local_seconds = f$local_seconds, score = apt_score(x, f$physical))
  }
  answer$cpu_seconds <- sum(proc.time()[1:2]) - cpu_before
  answer
}
aptopt_reset_contexts <- function() {
  # Reconstructible native pointers are never serialized into continuation state.
  assign(".aptopt_contexts", new.env(parent = emptyenv()), .GlobalEnv)
  invisible(NULL)
}
aptopt_pool <- function(x, workers, engine) {
  cl <- parallel::makePSOCKcluster(workers, outfile = "")
  ok <- FALSE; on.exit(if (!ok) parallel::stopCluster(cl))
  parallel::clusterCall(cl, function(wd) {
    setwd(wd); source("New_simulations/dm_apt_runtime_optimization/support.R"); NULL
  }, getwd())
  meta <- parallel::clusterCall(cl, aptopt_initialize, x, engine)
  apt_need(all(vapply(meta, function(v) identical(v$runtime, apt_runtime()), TRUE)), "Optimization worker runtime differs")
  ok <- TRUE; attr(cl, "runtime") <- meta; cl
}
aptopt_run_ladder <- function(x, plan, initial, engine, ladder = 1L, steps = plan$iterations,
                              restart_at = integer(), callback = NULL, pool = NULL) {
  s <- apt_new_state(x, plan, ladder, FALSE, initial)
  propagate <- function(requests) {
    requests <- lapply(requests, function(r) { r$engine <- engine; r })
    if (is.null(pool)) lapply(requests, aptopt_step) else parallel::parLapply(pool, requests, aptopt_step)
  }
  for (n in seq_len(steps)) {
    # EXACT original apt_iteration: one random adjacent exchange, then local
    # sweeps and unchanged diversity/temperature adaptation and bookkeeping.
    s <- apt_iteration(s, x, propagate)
    if (is.function(callback)) callback(n, s)
    if (n %in% restart_at) {
      s <- unserialize(serialize(s, NULL))
      if (is.null(pool)) aptopt_reset_contexts()
      else parallel::clusterCall(pool, aptopt_reset_contexts)
    }
  }
  s
}
