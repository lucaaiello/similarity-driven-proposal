# Deterministic fixtures and comparison helpers only. Sourcing this file does
# not run a sampler, optimizer, compiler, cluster worker, or checkpoint job.

aptopt_synthetic_fixture <- function(n = 16L, P = 8L, J = 4L) {
  stopifnot(n >= 8L, P >= 4L, J >= 3L)
  index <- seq_len(n * P)
  X <- unclass(scale(matrix(sin(index * .43) + cos(index * 1.37), n, P)))
  rownames(X) <- c("3026", "4014", as.character(seq.int(3L, n)))
  colnames(X) <- paste0("x", seq_len(P))
  B <- matrix(0, P, J); B[1, 1] <- .3; B[2, 2] <- -.25
  eta <- sweep(X %*% B, 2, seq(-.2, .2, length.out = J), "+")
  probability <- exp(eta) / rowSums(exp(eta))
  Y <- matrix(0L, n, J)
  for (i in seq_len(n)) {
    total <- 30L + i %% 9L
    counts <- floor(total * probability[i, ])
    remaining <- total - sum(counts)
    if (remaining > 0L) {
      order <- order(total * probability[i, ] - counts, decreasing = TRUE)
      counts[order[seq_len(remaining)]] <- counts[order[seq_len(remaining)]] + 1L
    }
    Y[i, ] <- as.integer(counts)
  }
  taxa <- c("Bacteroides", "Prevotella", paste0("Control", seq_len(J - 2L)))
  rates <- c(0, .05, .10, .15)
  start_names <- c("null", "moderate_5pct", "moderate_10pct", "dense_15pct")
  origins <- list()
  for (method in c("online_similarity_flip", "online_similarity_local", "uniform_flip")) {
    for (chain in seq_len(4L)) {
      xi <- matrix(0L, P, J)
      k <- as.integer(round(P * J * rates[chain]))
      if (k > 0L) {
        ordering <- order(sin(seq_len(P * J) * .79 + chain))
        xi[ordering[seq_len(k)]] <- 1L
      }
      state <- list(alpha = rep(0, J), xi = xi,
        beta = lapply(colSums(xi), function(k) numeric(k)),
        g = rep(1, J), Sigma_alpha = diag(J))
      origins[[sprintf("%s_chain_%02d", method, chain)]] <- list(method = method,
        chain = chain, origin_start = start_names[chain], start = state,
        origin_iteration = 0L, recorded_log_posterior = NULL)
    }
  }
  list(X = X, Y = Y, neighbors = lapply(seq_len(P), function(j) as.integer(setdiff(seq_len(P), j))),
    taxa = taxa, predictor_codes = colnames(X), start_names = start_names, origins = origins)
}

aptopt_worker_comparison <- function(reference, candidate, tolerance = 1e-10) {
  # Wall/CPU measurements are intentionally excluded from numerical parity.
  fields <- c("physical", "ll", "lp", "rng", "effective", "counts", "score")
  missing <- setdiff(fields, intersect(names(reference), names(candidate)))
  if (length(missing)) return(list(passed = FALSE, detail = paste("Missing fields:", paste(missing, collapse = ","))))
  exact <- list(rng = identical(reference$rng, candidate$rng),
    indicators = identical(reference$physical$xi, candidate$physical$xi),
    counts = identical(reference$counts, candidate$counts))
  numeric_fields <- c("ll", "lp", "effective", "score")
  numeric <- lapply(numeric_fields, function(field) all.equal(reference[[field]], candidate[[field]],
    tolerance = tolerance, check.attributes = TRUE))
  names(numeric) <- numeric_fields
  numeric$alpha <- all.equal(reference$physical$alpha, candidate$physical$alpha,
    tolerance = tolerance, check.attributes = TRUE)
  numeric$beta <- all.equal(reference$physical$beta, candidate$physical$beta,
    tolerance = tolerance, check.attributes = TRUE)
  list(passed = all(vapply(exact, isTRUE, TRUE)) && all(vapply(numeric, isTRUE, TRUE)),
    exact = exact, numeric = numeric)
}

aptopt_ladder_comparison <- function(reference, candidate, tolerance = 1e-10) {
  # Compare mathematical/Markov/adaptive state, never implementation receipts,
  # source hashes, compile times, process ids or serialized native pointers.
  fields <- c("completed", "plan", "slots", "rho", "master_rng", "labels", "seed_origin",
    "flow", "visits", "swaps_attempted", "swaps_accepted", "rho_projection",
    "beta_trace", "ll_trace", "lp_trace", "hypothetical_swap", "walker_trace",
    "scalar_trace", "g_trace", "counts", "local_crossings", "cold_swap_crossings",
    "roundtrip_trace", "cold_physical", "cold_fitted")
  missing <- setdiff(fields, intersect(names(reference), names(candidate)))
  if (length(missing)) return(list(passed = FALSE, detail = paste("Missing fields:", paste(missing, collapse = ","))))
  comparisons <- lapply(fields, function(field) all.equal(reference[[field]], candidate[[field]],
    tolerance = tolerance, check.attributes = TRUE))
  names(comparisons) <- fields
  exact <- list(master_rng = identical(reference$master_rng, candidate$master_rng),
    labels = identical(reference$labels, candidate$labels),
    walker_trace = identical(reference$walker_trace, candidate$walker_trace),
    swaps_attempted = identical(reference$swaps_attempted, candidate$swaps_attempted),
    swaps_accepted = identical(reference$swaps_accepted, candidate$swaps_accepted))
  for (k in seq_along(reference$slots)) {
    exact[[paste0("slot_rng_", k)]] <- identical(reference$slots[[k]]$rng, candidate$slots[[k]]$rng)
    exact[[paste0("slot_xi_", k)]] <- identical(reference$slots[[k]]$physical$xi, candidate$slots[[k]]$physical$xi)
  }
  list(passed = all(vapply(comparisons, isTRUE, TRUE)) && all(vapply(exact, isTRUE, TRUE)),
    fields = comparisons, exact = exact)
}
