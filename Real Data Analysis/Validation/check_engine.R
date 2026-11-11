# User-authorized bounded computational regression, never a production run.
# Modes parse/compile/fixed do not run MCMC;quick/full run explicitly bounded
# synthetic chains. Full additionally probes bounded real-data states/timings.
source("Real Data Analysis/Sampler/runtime.R")
apt_environment()
o <- apt_options(c("mode", "output", "workers"))
mode <- apt_get(o, "mode", "quick")
apt_need(mode %in% c("parse", "compile", "fixed", "quick", "full"), "Unknown optimization audit mode")
workers <- as.integer(apt_get(o, "workers", "2"))
apt_need(workers %in% c(1L, 2L, 5L, 10L), "Audit workers must be1,2,5or10")
for (f in aptopt_files()[grepl("[.]R$", aptopt_files())]) invisible(parse(f))
if (mode == "parse") { cat("OPTIMIZATION SOURCE PARSE PASSED; no compiler/optimizer/MCMC\n"); quit(save = "no") }
out <- apt_get(o,"output","Real Data Analysis/Validation/Results/Engine")
apt_need(!file.exists(out) && dir.create(out, recursive = TRUE), "Audit output already exists")
writeLines("INCOMPLETE", file.path(out, "status.txt"))
protected <- aptopt_protected(); sources <- aptopt_hashes(); checks <- list(); timings <- list()
check <- function(label, ok, detail = NULL) {
  checks[[length(checks) + 1L]] <<- list(check = label, passed = isTRUE(ok), detail = detail)
  if (!isTRUE(ok)) {
    saveRDS(checks, file.path(out, "checks.rds")); writeLines("FAILED", file.path(out, "status.txt"))
    stop("Persistent APT audit failed: ", label, "\n", paste(capture.output(str(detail)), collapse = "\n"), call. = FALSE)
  }
}
g <- readRDS("Real Data Analysis/Provenance/native_engine_receipt.rds")
check("generated_source_receipt", recorded_file_matches(g$reference_md5,"Real Data Analysis/Validation/Reference/Tempering/dm_sampler_reference.cpp") &&
  recorded_file_matches(g$native_md5,"Real Data Analysis/Sampler/dm_sampler.cpp"))
for (f in c("dm_sampler_reference.cpp", "joint_kernel.h", "dm_stable_numerics.h", "support.R", "math.R", "runner.R", "reference.R"))
  check(paste("reference_snapshot", f), identical(unname(tools::md5sum(file.path(apt_root, f))),
    unname(tools::md5sum(file.path("Real Data Analysis/Validation/Reference/Engine_Baseline", f)))))
x <- aptopt_synthetic_fixture(); aptopt_initialize(x, "both")
if (mode == "compile") {
  check("protected_sources_unchanged", identical(protected, aptopt_protected()))
  saveRDS(checks, file.path(out, "checks.rds")); writeLines("COMPILE_ONLY", file.path(out, "status.txt"))
  cat("BOTH OPTIMIZATION ENGINES COMPILED;no optimizer/MCMC\n"); quit(save = "no")
}
reference <- get(".apt_cpp", .GlobalEnv); candidate <- get(".aptopt_cpp", .GlobalEnv)
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(2026102901L)
seed <- .Random.seed
for (method in c("online_similarity_local", "online_similarity_flip", "uniform_flip")) {
  ctx <- aptopt_context(method); info <- candidate$apt_context_info(ctx)
  check(paste("constructor_no_rng", method), identical(seed, .Random.seed))
  check(paste("constructor_data_dimensions", method), info$n == nrow(x$X) && info$P == ncol(x$X) && info$J == ncol(x$Y))
  broken <- unserialize(serialize(ctx, NULL))
  check(paste("serialized_pointer_rejected", method), inherits(try(candidate$apt_context_info(broken), silent = TRUE), "try-error"))
}
# Prefix target/proposal helpers remain numerically unchanged. No RNG or fits.
xi <- matrix(0L, ncol(x$X), ncol(x$Y)); xi[1:2, 2] <- 1L
B <- matrix(0, ncol(x$X), ncol(x$Y)); B[1:2, 2] <- c(.1, -.1)
alpha <- rep(.1, ncol(x$Y)); cn <- c(.12, .11, -.08)
for (t in c(1, .4, .01)) {
  args <- list(x$X[, 1:2], x$Y, alpha, x$X %*% B, 2L, B[1:2, 2], cn, 10, 10, 2, t)
  a <- do.call(reference$audit_dm_joint_block_terms, args)
  b <- do.call(candidate$audit_dm_joint_block_terms, args)
  check(paste("fixed_joint_target_precision", t), isTRUE(all.equal(a, b, tolerance = 1e-12)))
}
if (mode %in% c("quick", "full")) {
  cat("BOUNDED COMPUTATIONAL TEST:synthetic MCMC will run;no scientific production chain or cluster job.\n")
  stream <- .Random.seed
  for (method in c("online_similarity_local", "online_similarity_flip", "uniform_flip"))
    for (origin in 1:4) for (t in c(1, .4, .01)) {
      start <- x$origins[[sprintf("%s_chain_%02d", method, origin)]]$start
      stream <- parallel::nextRNGStream(stream)
      request <- list(physical = apt_physical(start), tuning = start[c("g", "Sigma_alpha")],
        rng = stream, beta = t, method = method)
      request$engine <- "reference"; a <- aptopt_step(request)
      request$engine <- "persistent"; b <- aptopt_step(request)
      result <- aptopt_worker_comparison(a, b)
      check(paste("fresh_worker", method, origin, t), result$passed, result)
    }
  z <- apt_load(); p <- z$plan
  p$replicas <- if (mode == "full") 10L else 3L
  p$iterations <- if (mode == "full") 80L else 24L
  p$burn <- if (mode == "full") 20L else 8L
  p$thin <- 2L; p$adapt_until <- p$iterations; p$g_window <- 5L
  initial <- apt_worker_initial_scores(p$method)
  reference_states <- vector("list", p$iterations)
  a <- aptopt_run_ladder(x, p, initial, "reference", callback = function(n, q) reference_states[[n]] <<- q)
  aptopt_reset_contexts()
  b <- aptopt_run_ladder(x, p, initial, "persistent", restart_at = c(7L, 17L), callback = function(n, q) {
    r <- aptopt_ladder_comparison(reference_states[[n]], q)
    check(paste("whole_ladder_every_boundary", n), r$passed, r)
  })
  check("whole_ladder_has_accepted_swaps", sum(a$swaps_accepted) > 0)
  check("whole_ladder_g_updates_occur", all(vapply(a$slots, function(q) q$tuning$g_epoch > 0, TRUE)))
  r <- aptopt_ladder_comparison(a, b); check("whole_ladder_serialized_restart", r$passed, r)
  # Scheduling and different worker counts may not alter slot RNG or decisions.
  cl <- aptopt_pool(x, workers, "persistent")
  parallel_result <- tryCatch(aptopt_run_ladder(x, p, initial, "persistent", pool = cl,
    restart_at = c(7L, 17L)), finally = parallel::stopCluster(cl))
  r <- aptopt_ladder_comparison(a, parallel_result)
  check(paste("parallel_whole_ladder_workers", workers), r$passed, r)
  rm(reference_states, a, b, parallel_result); gc(FALSE)
}
if (mode == "full") {
  # Explicitly bounded full-size one-step comparisons. Saved states are diagnostic
  # fixtures only; they are never a production initialization recommendation.
  real <- apt_load()$inputs
  assign(".apt_x", real, .GlobalEnv); aptopt_reset_contexts()
  plan <- apt_load()$plan; real_cases <- list()
  for (method in c("online_similarity_local", "online_similarity_flip", "uniform_flip")) {
    for (chain in c(1L, 3L)) {
      start <- real$origins[[sprintf("%s_chain_%02d", method, chain)]]$start
      real_cases[[length(real_cases) + 1L]] <- list(physical = apt_physical(start),
        tuning = start[c("g", "Sigma_alpha")], beta = if (chain == 1L) 1 else .3, method = method)
    }
  }
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(2026102902L); stream <- .Random.seed
  for (i in seq_along(real_cases)) {
    request <- real_cases[[i]]; stream <- parallel::nextRNGStream(stream); request$rng <- stream
    order <- if (i %% 2) c("reference", "persistent") else c("persistent", "reference")
    results <- list()
    for (engine in order) {
      request$engine <- engine; begin <- proc.time()[["elapsed"]]
      results[[engine]] <- aptopt_step(request)
      timings[[length(timings) + 1L]] <- data.frame(case = i, method = request$method,
        engine = engine, seconds = proc.time()[["elapsed"]] - begin,
        native_seconds = results[[engine]]$local_seconds)
    }
    r <- aptopt_worker_comparison(results$reference, results$persistent)
    check(paste("full_size_saved_worker", i), r$passed, r)
  }
}
check("protected_original_sources_results_unchanged", identical(protected, aptopt_protected()))
check("optimization_sources_unchanged_during_test", identical(sources, aptopt_hashes()))
saveRDS(checks, file.path(out, "checks.rds"))
timing <- if (length(timings)) do.call(rbind, timings) else data.frame()
saveRDS(list(passed = TRUE, mode = mode, workers = workers, checks = length(checks),
  source_hashes = sources, original_hashes = protected, runtime = apt_runtime(),
  timing = timing, schedule = "unchanged_one_uniformly_random_adjacent_exchange_per_iteration"),
  file.path(out, "audit.rds"))
report <- c(paste("PERSISTENT APT", toupper(mode), "VERIFICATION PASSED"), paste("Checks:", length(checks)),
  "Original one-random-adjacent-swap schedule,23.4% target,temperature/diversity adaptation and local kernels preserved.",
  if (mode %in% c("quick", "full")) paste("Bounded synthetic whole-ladder tests require exact indicators,RNG,swap decisions and numeric agreement1e-10;workers", workers)
  else "No sampler or optimizer executed.",
  if (mode == "full") "Six additional real-data single-step FIXTURES compared;not production runs or converged inference.",
  "Original sampler,source graph,plans,completed checkpoints and results unchanged.",
  "Timing probes exclude compilation,rotate engine order,and do not establish a whole-chain or cluster speedup.",
  if (nrow(timing)) capture.output(print(timing, row.names = FALSE, digits = 5)))
writeLines(report, file.path(out, "report.txt")); writeLines("COMPLETE", file.path(out, "status.txt"))
cat(report, sep = "\n"); cat("\nOutput:", out, "\n")
