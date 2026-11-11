# Orchestration only. Each PSOCK worker has one OpenMP/BLAS thread and a private
# compiled-module cache. Complete physical states swap; proposal tuning and RNG
# streams stay at temperature slots. All state needed for restart is serialized.
apt_new_state <- function(x, plan, ladder, control, initial) {
  L <- if (control) 1L else plan$replicas
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(plan$seed_base + ladder)
  master_rng <- .Random.seed; stream <- master_rng
  seed_origin <- if (control) ladder else (ladder - 1L + seq_len(L) - 1L) %% 4L + 1L
  origin <- x$origins[[sprintf("%s_chain_%02d", plan$method, ladder)]]$start
  slots <- lapply(seq_len(L), function(k) {
    stream <<- parallel::nextRNGStream(stream)
    c(initial[[seed_origin[k]]], list(rng = stream,
      tuning = list(g = origin$g, Sigma_alpha = origin$Sigma_alpha,
        effective_sum = rep(0, ncol(x$Y)), effective_count = 0L, g_epoch = 0L)))
  })
  N <- plan$iterations; nsave <- floor((N - plan$burn) / plan$thin)
  edges <- max(0L, L - 1L)
  list(schema = "dm_mmv_pt_checkpoint_v1", completed = 0L, ladder = ladder, control = control,
    plan = plan, source_hashes = apt_hashes(),
    input_hash = tools::md5sum("Real Data Analysis/Data/Reference/reference_plan.rds"), runtime = apt_runtime(),
    slots = slots, rho = if (L > 1L) rep(log(log(plan$initial_tmax) / (L - 1L)), L - 1L) else numeric(),
    master_rng = master_rng, labels = seq_len(L), seed_origin = seed_origin,
    flow = list(phase = replace(integer(L), 1L, 1L), departed = integer(L),
      roundtrips = integer(L), durations = integer()), visits = matrix(0L, L, L),
    swaps_attempted = matrix(0L, 2L, edges), swaps_accepted = matrix(0L, 2L, edges),
    rho_projection = integer(edges),
    beta_trace = matrix(NA_real_, N, L), ll_trace = matrix(NA_real_, N, L),
    lp_trace = matrix(NA_real_, N, L), hypothetical_swap = matrix(NA_real_, N, edges),
    walker_trace = matrix(NA_integer_, N, L),
    scalar_trace = array(NA_real_, c(N, L, 6L), dimnames = list(NULL, NULL,
      c("prev_alpha", "prev_size", "model_size", "subject3026", "subject4014", "mean_log_total"))),
    g_trace = array(NA_real_, c(N, L, ncol(x$Y))),
    counts = matrix(0, L, 5L, dimnames = list(NULL, c("alpha", "flip", "local", "joint", "projections"))),
    native_seconds = numeric(L), worker_cpu_seconds = numeric(L),
    local_crossings = array(0L, c(L, 2L, 2L)), cold_swap_crossings = matrix(0L, 2L, 2L),
    roundtrip_trace = integer(N),
    cold_physical = vector("list", nsave), cold_fitted = matrix(NA_real_, nrow(x$Y), nsave),
    elapsed_seconds = 0, startup_seconds = 0, allocated_slot_seconds = 0)
}
apt_swap_slots <- function(s, pair) {
  for (name in c("physical", "ll", "lp", "score")) {
    temp <- s$slots[[pair]][[name]]
    s$slots[[pair]][[name]] <- s$slots[[pair + 1L]][[name]]
    s$slots[[pair + 1L]][[name]] <- temp
  }
  s$labels[c(pair, pair + 1L)] <- s$labels[c(pair + 1L, pair)]
  s
}
apt_iteration <- function(s, x, propagate) {
  n <- s$completed + 1L; L <- length(s$slots); p <- s$plan
  phase <- if (n <= p$burn) 1L else 2L
  beta <- if (L == 1L) 1 else apt_beta(s$rho)
  id <- match(c("3026", "4014"), rownames(x$X))
  assign(".Random.seed", s$master_rng, envir = .GlobalEnv)
  # MMV ordering: one uniformly chosen adjacent swap, then local propagation,
  # then ALL adjacent-pair probability signals for temperature adaptation.
  if (L > 1L) {
    prob <- apt_swap_prob(vapply(s$slots, function(z) z$ll, 0), beta)
    pair <- sample.int(L - 1L, 1L)
    s$swaps_attempted[phase, pair] <- s$swaps_attempted[phase, pair] + 1L
    old_cold <- s$slots[[1]]$score$fitted[id] > .05
    if (runif(1) < prob[pair]) {
      s <- apt_swap_slots(s, pair)
      s$swaps_accepted[phase, pair] <- s$swaps_accepted[phase, pair] + 1L
      s$cold_swap_crossings[phase, ] <- s$cold_swap_crossings[phase, ] +
        as.integer((s$slots[[1]]$score$fitted[id] > .05) != old_cold)
    }
  }
  s$master_rng <- get(".Random.seed", envir = .GlobalEnv)
  requests <- lapply(seq_len(L), function(k) list(physical = s$slots[[k]]$physical,
    tuning = s$slots[[k]]$tuning, rng = s$slots[[k]]$rng, beta = beta[k], method = p$method))
  values <- propagate(requests)
  for (k in seq_len(L)) {
    old <- s$slots[[k]]; v <- values[[k]]
    s$local_crossings[k, , phase] <- s$local_crossings[k, , phase] +
      as.integer((old$score$fitted[id] > .05) != (v$score$fitted[id] > .05))
    tuning <- apt_g_update(old$tuning, v$effective, p, p$method != "uniform_flip")
    s$slots[[k]] <- c(v[c("physical", "ll", "lp", "rng", "score")], list(tuning = tuning))
    s$counts[k, ] <- s$counts[k, ] + v$counts
    s$native_seconds[k] <- s$native_seconds[k] + v$local_seconds
    s$worker_cpu_seconds[k] <- s$worker_cpu_seconds[k] + v$cpu_seconds
    s$ll_trace[n, k] <- v$ll; s$lp_trace[n, k] <- v$lp
    s$scalar_trace[n, k, ] <- c(v$score$prev_alpha, v$score$prev_size, v$score$model_size,
      v$score$fitted[id], v$score$mean_log_total)
    s$g_trace[n, k, ] <- tuning$g
  }
  if (L > 1L) {
    update <- apt_rho_update(s$rho, vapply(s$slots, function(v) v$ll, 0), n,
      p$target_swap, p$gain, p$decay, p$rho_lower, p$rho_upper)
    s$hypothetical_swap[n, ] <- update$probabilities
    if (n <= p$adapt_until) {
      s$rho <- update$rho; s$rho_projection <- s$rho_projection + as.integer(update$projected)
    }
    s$flow <- apt_flow_update(s$flow, s$labels, n)
  }
  s$beta_trace[n, ] <- beta; s$walker_trace[n, ] <- s$labels
  s$roundtrip_trace[n] <- sum(s$flow$roundtrips)
  for (k in seq_len(L)) s$visits[s$labels[k], k] <- s$visits[s$labels[k], k] + 1L
  if (n > p$burn && (n - p$burn) %% p$thin == 0L) {
    i <- (n - p$burn) %/% p$thin
    s$cold_physical[[i]] <- s$slots[[1]]$physical
    s$cold_fitted[, i] <- s$slots[[1]]$score$fitted
  }
  s$completed <- n; s
}
apt_run <- function(x, plan, ladder, control, out, workers, audit = FALSE) {
  apt_need(ladder %in% 1:4, "Ladder must be1:4")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(out, "task_lock")
  apt_need(dir.create(lock, showWarnings = FALSE), "Task lock exists; verify its job is dead before recovery")
  writeLines(c(paste("job", Sys.getenv("JOB_ID", "local")), paste("pid", Sys.getpid())), file.path(lock, "owner.txt"))
  cl <- NULL
  on.exit({
    if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
    file.remove(file.path(lock, "owner.txt"))
    if (identical(basename(lock), "task_lock") &&
      identical(normalizePath(dirname(lock)), normalizePath(out)) &&
      !length(list.files(lock, all.files = TRUE, no.. = TRUE))) unlink(lock, recursive = TRUE)
  })
  t0 <- proc.time()[["elapsed"]]
  checkpoints <- list.dirs(file.path(out, "checkpoints"), full.names = TRUE, recursive = FALSE)
  checkpoints <- checkpoints[grepl("^iter_[0-9]{6}$", basename(checkpoints))]
  if (length(checkpoints)) {
    s <- apt_read(tail(sort(checkpoints), 1))
    apt_need(identical(s$plan, plan) && identical(s$source_hashes, apt_hashes()) &&
      identical(s$input_hash, tools::md5sum("Real Data Analysis/Data/Reference/reference_plan.rds")) &&
      identical(s$runtime, apt_runtime()) && s$ladder == ladder && identical(s$control, control),
      "Cannot resume: checkpoint provenance/plan/runtime differs")
    cat("RESUMING complete checkpoint at", s$completed, "\n")
  } else s <- NULL
  if (!is.null(s) && s$completed == plan$iterations) { cat("KEEPING completed task\n"); return(invisible(s)) }
  cl <- apt_pool(x, workers)
  if (is.null(s)) {
    initial <- parallel::clusterCall(cl, apt_worker_initial_scores, plan$method)[[1]]
    s <- apt_new_state(x, plan, ladder, control, initial)
  }
  startup <- proc.time()[["elapsed"]] - t0
  prior_elapsed <- s$elapsed_seconds; prior_slots <- s$allocated_slot_seconds
  s$startup_seconds <- s$startup_seconds + startup
  cat("START", if (control) "UNTEMPERED CONTROL" else "ADAPTIVE TEMPERING", "ladder", ladder,
    "replicas", length(s$slots), "workers", workers, "from", s$completed, "to", plan$iterations, "\n")
  propagate <- function(requests) parallel::parLapply(cl, requests, apt_worker_step)
  segment_start <- s$completed
  for (n in seq.int(s$completed + 1L, plan$iterations)) {
    s <- apt_iteration(s, x, propagate)
    s$elapsed_seconds <- prior_elapsed + proc.time()[["elapsed"]] - t0
    s$allocated_slot_seconds <- prior_slots + workers * (proc.time()[["elapsed"]] - t0)
    if (n %% 50L == 0L || n == 1L) {
      current_work <- max(.001, proc.time()[["elapsed"]] - t0 - startup)
      # Segment-aware ETA after restart; never divide by the old completed count.
      segment_done <- n - segment_start
      eta <- current_work / segment_done * (plan$iterations - n)
      cat(sprintf("APT iteration=%d/%d beta_hot=%.6g roundtrips=%d elapsed=%.1fs eta=%.1fs\n",
        n, plan$iterations, if (length(s$rho)) tail(apt_beta(s$rho), 1) else 1,
        sum(s$flow$roundtrips), s$elapsed_seconds, eta)); flush.console()
    }
    if (n %% plan$checkpoint_every == 0L || n == plan$iterations) {
      apt_need(identical(s$source_hashes, apt_hashes()), "Sources changed during APT task")
      apt_write(s, file.path(out, "checkpoints", sprintf("iter_%06d", n)))
      cat("COMMITTED APT CHECKPOINT", n, "\n"); flush.console()
    }
  }
  writeLines("COMPLETE", file.path(out, "status.txt"))
  invisible(s)
}
