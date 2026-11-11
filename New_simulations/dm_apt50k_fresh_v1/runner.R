# Existing local/adaptation bookkeeping; explicit new exchange schedule.
long_iteration <- function(s, x, propagate) {
  n <- s$completed + 1L; L <- length(s$slots); p <- s$plan
  phase <- if (n <= p$burn) 1L else 2L
  beta <- if (L == 1L) 1 else apt_beta(s$rho)
  id <- match(c("3026", "4014"), rownames(x$X))
  assign(".Random.seed", s$master_rng, envir = .GlobalEnv)
  # One random-order sweep of all adjacent exchange kernels, then local
  # propagation and the existing all-pair temperature-adaptation signals.
  if (L > 1L) s <- long_exchange(s, x, n, beta)
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

long_trim_previous <- function(previous, out, milestones) {
  if (is.null(previous) || !dir.exists(previous)) return(invisible(NULL))
  iteration <- as.integer(sub("iter_", "", basename(previous), fixed = TRUE))
  if (iteration %in% milestones) return(invisible(NULL))
  # Only this NEW run's completed, superseded rolling checkpoint may be removed.
  parent <- normalizePath(file.path(out, "checkpoints"), winslash = "/", mustWork = TRUE)
  candidate <- normalizePath(previous, winslash = "/", mustWork = TRUE)
  apt_need(identical(dirname(candidate), parent) && grepl("^iter_[0-9]{6}$", basename(candidate)) &&
    identical(sort(list.files(candidate, all.files = TRUE, no.. = TRUE)),
      sort(c("state.rds", "state.md5", "status.txt"))) &&
    identical(readLines(file.path(candidate, "status.txt")), "COMPLETE"), "Unsafe rolling-checkpoint cleanup refused")
  unlink(candidate, recursive = TRUE)
  apt_need(!dir.exists(candidate), "Cannot remove superseded rolling checkpoint")
}
long_run <- function(x, plan, start, out, workers) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(out, "task_lock")
  apt_need(dir.create(lock), "Task lock exists; inspect before retrying")
  writeLines(c(paste("pid", Sys.getpid()), paste("job", Sys.getenv("JOB_ID", "local"))), file.path(lock, "owner.txt"))
  cl <- NULL
  on.exit({
    if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
    file.remove(file.path(lock, "owner.txt"))
    if (!length(list.files(lock, all.files = TRUE, no.. = TRUE))) unlink(lock, recursive = TRUE)
  })
  t0 <- proc.time()[["elapsed"]]
  checkpoints <- list.dirs(file.path(out, "checkpoints"), recursive = FALSE, full.names = TRUE)
  checkpoints <- sort(checkpoints[grepl("^iter_[0-9]{6}$", basename(checkpoints))])
  previous <- if (length(checkpoints)) tail(checkpoints, 1) else NULL
  s <- if (is.null(previous)) NULL else apt_read(previous)
  if (!is.null(s)) {
    apt_need(identical(s$schema, "dm_apt50k_fresh_v1") && identical(s$plan, plan) &&
      identical(s$source_hashes, long_hashes()) && identical(s$input_hash,
        tools::md5sum(file.path(long_root, "inputs/plan.rds"))) && identical(s$runtime, apt_runtime()) &&
      s$ladder == start, "Resume plan/source/input/runtime differs")
    if (s$completed == plan$iterations) { cat("Keeping completed task\n"); return(invisible(s)) }
  }
  cl <- aptopt_pool(x, workers, "persistent")
  if (is.null(s)) {
    initial <- parallel::clusterCall(cl, function(method, start) {
      source("New_simulations/dm_apt50k_fresh_v1/support.R")
      long_seed_scores(get(".apt_x", .GlobalEnv), method, start)
    }, plan$method, start)[[1]]
    s <- long_new(x, plan, start, initial)
  }
  startup <- proc.time()[["elapsed"]] - t0
  old_elapsed <- s$elapsed_seconds; old_allocated <- s$allocated_slot_seconds
  s$startup_seconds <- s$startup_seconds + startup; from <- s$completed
  propagate <- function(requests) parallel::parLapply(cl, requests, aptopt_step)
  cat("START FRESH50K", plan$method, "start", long_names[start], "replicas10 workers10 from", from,
    "to", plan$iterations, "burn", plan$burn, "exchanges9/sweep\n"); flush.console()
  for (n in seq.int(from + 1L, plan$iterations)) {
    s <- long_iteration(s, x, propagate)
    now <- proc.time()[["elapsed"]] - t0
    s$elapsed_seconds <- old_elapsed + now
    s$allocated_slot_seconds <- old_allocated + workers * now
    if (n == from + 1L || n %% 50L == 0L) {
      eta <- max(.001, now - startup) / (n - from) * (plan$iterations - n)
      cat(sprintf("FRESH50K method=%s start=%d iteration=%d/%d beta_hot=%.6g roundtrips=%d elapsed=%.1fs eta=%.1fs\n",
        plan$method, start, n, plan$iterations, tail(apt_beta(s$rho), 1), sum(s$flow$roundtrips),
        s$elapsed_seconds, eta)); flush.console()
    }
    if (n %% plan$checkpoint_every == 0L || n == plan$iterations) {
      long_gate()
      dest <- file.path(out, "checkpoints", sprintf("iter_%06d", n))
      apt_write(s, dest)
      cat("COMMITTED FRESH50K CHECKPOINT", n, "\n"); flush.console()
      long_trim_previous(previous, out, plan$milestones); previous <- dest
    }
  }
  writeLines("COMPLETE", file.path(out, "status.txt"))
  invisible(s)
}
