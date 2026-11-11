args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

mode <- arg_value("mode", "quick")
stage <- arg_value("stage", "all")
skip_audit <- identical(tolower(arg_value("skip-audit", "false")), "true")
if (!mode %in% c("quick", "full")) stop("mode must be quick or full")
if (!stage %in% c("all", "core", "plugin")) stop("stage must be all, core, or plugin")

datasets <- as.integer(arg_value("datasets", if (mode == "quick") "1" else "20"))
chains <- as.integer(arg_value("chains", if (mode == "quick") "2" else "4"))
runs <- as.integer(arg_value("runs", if (mode == "quick") "300" else "20000"))
burn <- as.integer(arg_value("burn", if (mode == "quick") "150" else "10000"))
thin <- as.integer(arg_value("thin", if (mode == "quick") "1" else "5"))
bool_option <- function(name) {
  value <- tolower(arg_value(name,"false"))
  if (!value %in% c("true","false")) stop(name," must be true or false")
  identical(value,"true")
}
coefficient_warm_start <- bool_option("coefficient-warm-start")
within_model_update <- bool_option("within-model-update")
save_full_fits <- bool_option("save-full-fits")
score_maxit_core <- as.integer(arg_value("score-optimizer-maxit", "1"))
cache_option <- tolower(arg_value("cache-score-fits", "auto"))
if (!cache_option %in% c("auto", "true", "false") ||
    is.na(score_maxit_core) || score_maxit_core < 1L) stop("Invalid score fitting settings")
core_cache_override <- if (cache_option == "auto") NULL else cache_option == "true"
coefficient_maxit <- as.integer(arg_value("coefficient-optimizer-maxit","3"))
within_attempts <- as.integer(arg_value("within-model-attempts","1"))
within_scale <- as.numeric(arg_value("within-model-scale","0.25"))
if (anyNA(c(coefficient_maxit,within_attempts,within_scale)) || coefficient_maxit<1 ||
    within_attempts<1 || !is.finite(within_scale) || within_scale<=0) stop("Invalid coefficient settings")
if (coefficient_warm_start || within_model_update) {
  source("Simulations/Validation/dm_validation.R")
  gate_dir <- "Simulations/Validation/Results/New_Run/DM_Coefficients"
  gate <- read.csv(file.path(gate_dir,"summary.csv"))
  gate_hash <- read.csv(file.path(gate_dir,"source_md5.csv"))
  hashes <- tools::md5sum(dm_source_files)
  if (nrow(gate)!=1L || !isTRUE(gate$passed[1]) ||
      !identical(gate_hash$file,names(hashes)) || !identical(gate_hash$md5,unname(hashes)))
    stop("Optional coefficient-update audit missing/stale/failed")
  if (!length(grep("^--output=",args))) stop("Use an explicit fresh --output for optional-kernel experiments")
}
dataset_ids_text <- arg_value("dataset-ids", "")
method_names_text <- arg_value("methods", "")
out_dir <- arg_value(
  "output",
  "Simulations/Dirichlet_Multinomial/Comparison/Results/New_Run"
)
if (datasets < 1L || chains < 2L || runs <= burn || burn < 0L || thin < 1L) {
  stop("invalid datasets, chains, runs, burn, or thin")
}
dataset_ids <- if (nzchar(dataset_ids_text)) {
  as.integer(strsplit(dataset_ids_text, ",", fixed = TRUE)[[1]])
} else {
  seq_len(datasets)
}
if (!length(dataset_ids) || any(!is.finite(dataset_ids)) || any(dataset_ids < 1L)) {
  stop("dataset-ids must be a comma-separated list of positive integers")
}

if (mode == "full" && !skip_audit) {
  audit_file <- file.path("Simulations/Validation/Results/New_Run/DM_Adaptation/summary.csv")
  if (!file.exists(audit_file)) {
    stop("Run Simulations/Validation/check_dm_adaptation.R successfully before the full experiment")
  }
  audit <- read.csv(audit_file, stringsAsFactors = FALSE)
  if (!nrow(audit) || !isTRUE(as.logical(audit$passed[[1]]))) {
    stop("the DM online-adaptation audit has not passed")
  }
}

if ((save_full_fits || score_maxit_core != 1L || cache_option != "auto") &&
    dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Use a fresh --output directory for diagnostic/sensitivity runs")
}
if (mode == "full") {
  source("Simulations/Validation/dm_validation.R")
  stability_dir <- "Simulations/Validation/Results/New_Run/DM_Numerics"
  stability <- read.csv(file.path(stability_dir, "summary.csv"))
  stability_hash <- read.csv(file.path(stability_dir, "source_md5.csv"))
  hashes <- tools::md5sum(dm_source_files)
  if (nrow(stability) != 1L || !isTRUE(stability$passed[1]) ||
      !identical(stability_hash$file, names(hashes)) ||
      !identical(stability_hash$md5, unname(hashes))) {
    stop("Run Simulations/Validation/check_dm_numerics.R: current-source numerical audit missing/stale/failed")
  }
}
if (save_full_fits && mode != "full") source("Simulations/Validation/dm_validation.R")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (save_full_fits) dir.create(file.path(out_dir, "full_fits"), showWarnings = FALSE)
source_files <- c("Simulations/Dirichlet_Multinomial/run_dirichlet_multinomial.R",
                  "Simulations/Dirichlet_Multinomial/dm_sampler.cpp",
                  "Simulations/Dirichlet_Multinomial/dm_numerics.h", "Simulations/Validation/dm_validation.R")
source_hashes <- tools::md5sum(source_files)
write.csv(data.frame(file = names(source_hashes), md5 = unname(source_hashes)),
          file.path(out_dir, "source_md5.csv"), row.names = FALSE)
snapshot_dir <- file.path(out_dir, "source_snapshot")
dir.create(snapshot_dir, showWarnings = FALSE)
if (!all(file.copy(source_files, snapshot_dir, overwrite = FALSE))) {
  stop("Could not archive generating source; use a fresh --output directory")
}
cache_dir <- file.path("Simulations/.rcpp_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(
  file.path("Simulations/Dirichlet_Multinomial/dm_sampler.cpp"),
  cacheDir = cache_dir,
  rebuild = FALSE,
  showOutput = FALSE,
  verbose = FALSE
)

ess_basic <- function(x) {
  if (length(x) < 20L || !is.finite(var(x)) || var(x) == 0) return(NA_real_)
  z <- as.numeric(acf(x, lag.max = min(length(x) - 1L, 1000L), plot = FALSE)$acf)[-1]
  m <- length(z) - length(z) %% 2L
  if (m < 2L) return(length(x))
  pairs <- z[seq(1L, m, 2L)] + z[seq(2L, m, 2L)]
  first_negative <- match(TRUE, pairs <= 0, nomatch = length(pairs) + 1L)
  keep <- if (first_negative == 1L) numeric() else pairs[seq_len(first_negative - 1L)]
  min(length(x), length(x) / (1 + 2 * sum(keep)))
}

bulk_ess <- function(x) {
  if (length(x) < 20L) return(NA_real_)
  z <- qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
  ess_basic(z)
}

tail_ess <- function(x) {
  if (length(x) < 20L) return(NA_real_)
  limits <- quantile(x, c(0.05, 0.95), names = FALSE)
  values <- c(ess_basic(as.numeric(x <= limits[1])),
              ess_basic(as.numeric(x >= limits[2])))
  if (all(is.na(values))) NA_real_ else min(values, na.rm = TRUE)
}

rank_split_rhat <- function(chain_list) {
  common_length <- min(vapply(chain_list, length, integer(1)))
  half <- floor(common_length / 2)
  if (half < 20L) return(NA_real_)
  split_chains <- unlist(lapply(chain_list, function(x) {
    list(x[seq_len(half)], x[(common_length - half + 1L):common_length])
  }), recursive = FALSE)
  pooled <- unlist(split_chains, use.names = FALSE)
  z <- qnorm((rank(pooled, ties.method = "average") - 0.5) / length(pooled))
  draws <- matrix(z, nrow = half)
  within <- mean(apply(draws, 2L, var))
  if (!is.finite(within) || within == 0) return(1)
  between <- half * var(colMeans(draws))
  sqrt((((half - 1) / half) * within + between / half) / within)
}

simulate_dm <- function(seed, n = 70L, J = 5L, P = 18L, rho = 0.65,
                        total_count = 100L) {
  set.seed(seed)
  sigma <- outer(seq_len(P), seq_len(P), function(i, j) rho^abs(i - j))
  X <- matrix(rnorm(n * P), n, P) %*% chol(sigma)
  X <- scale(X)
  beta <- matrix(0, P, J)
  active <- rbind(
    c(1, 1), c(4, 1), c(2, 2), c(7, 2), c(5, 3),
    c(10, 3), c(8, 4), c(13, 4), c(11, 5), c(16, 5)
  )
  effects <- c(0.70, -0.60, -0.65, 0.55, 0.65, -0.55, -0.70, 0.60, 0.55, -0.60)
  beta[active] <- effects
  alpha <- log(c(1.4, 1.1, 0.9, 1.2, 0.8))
  eta <- sweep(X %*% beta, 2L, alpha, "+")
  # Scalar-first pmin/pmax can drop matrix dimensions.  Clamp in place so the
  # observation-by-category structure is preserved for row-wise DM sampling.
  eta[eta < -5] <- -5
  eta[eta > 5] <- 5
  gamma <- exp(eta)
  Y <- matrix(0L, n, J)
  for (i in seq_len(n)) {
    latent <- rgamma(J, shape = gamma[i, ], rate = 1)
    probabilities <- latent / sum(latent)
    Y[i, ] <- as.integer(rmultinom(1L, total_count, probabilities))
  }
  list(X = unclass(X), Y = Y, alpha = alpha, beta = beta, xi = (beta != 0))
}

make_neighbors <- function(X, degree = 4L) {
  correlation <- abs(cor(X))
  diag(correlation) <- -Inf
  lapply(seq_len(ncol(X)), function(j) {
    order(correlation[j, ], decreasing = TRUE)[seq_len(min(degree, ncol(X) - 1L))]
  })
}

make_start <- function(chain, truth) {
  P <- nrow(truth)
  J <- ncol(truth)
  if (chain == 1L) return(matrix(0L, P, J))
  if (chain == 2L) return(matrix(rbinom(P * J, 1L, 0.05), P, J))
  if (chain == 3L) return(matrix(rbinom(P * J, 1L, 0.20), P, J))
  start <- truth * 1L
  switch <- sample.int(P * J, max(2L, floor(0.10 * P * J)))
  start[switch] <- 1L - start[switch]
  start
}

full_beta_draws <- function(fit) {
  dims <- dim(fit$xi)
  result <- array(0, dims)
  for (draw in seq_len(dims[3])) {
    for (category in seq_len(dims[2])) {
      selected <- which(fit$xi[, category, draw] == 1L)
      if (length(selected)) result[selected, category, draw] <- fit$beta[[draw]][[category]]
    }
  }
  result
}

all_core_method_names <- c(
  "uniform_flip", "legacy_similarity_flip", "optimized_similarity_fixed",
  "online_similarity_flip", "online_similarity_local", "gagnon_laplace_rj"
)
core_specs <- data.frame(
  method = all_core_method_names,
  informed = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  adapt = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
  balanced = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
  mixture = c(0, 0, 0.02, 0.02, 0.02, 0),
  local_move = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
  cache = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  proposal = c(rep("similarity", 5L), "gagnon_laplace"),
  laplace_maxit = c(rep(25L, 5L), 25L),
  stringsAsFactors = FALSE
)
if (nzchar(method_names_text)) {
  method_names <- strsplit(method_names_text, ",", fixed = TRUE)[[1]]
  unknown_methods <- setdiff(method_names, core_specs$method)
  if (length(unknown_methods)) {
    stop("unknown methods: ", paste(unknown_methods, collapse = ", "))
  }
  core_specs <- core_specs[core_specs$method %in% method_names, , drop = FALSE]
}
if (mode == "full" && !skip_audit && "gagnon_laplace_rj" %in% core_specs$method) {
  gagnon_audit_file <- file.path(
    "Simulations/Validation/Results/New_Run/DM_Laplace/summary.csv"
  )
  if (!file.exists(gagnon_audit_file)) {
    stop("Run Simulations/Validation/check_dm_laplace.R successfully before the full Gagnon comparison")
  }
  gagnon_audit <- read.csv(gagnon_audit_file, stringsAsFactors = FALSE)
  if (!nrow(gagnon_audit) || !isTRUE(as.logical(gagnon_audit$passed[[1]]))) {
    stop("the DM Gagnon-Laplace proposal audit has not passed")
  }
}

chain_rows <- list()
pip_rows <- list()
coefficient_rows <- list()
trace_rows <- list()
adaptation_rows <- list()
failures <- list()
chain_id <- pip_id <- coefficient_id <- trace_id <- adaptation_id <- failure_id <- 0L

run_one <- function(dat, dataset, chain, spec, score_maxit = 1L,
                    quality = "production", cache_override = NULL) {
  method_seed_index <- match(spec$method, all_core_method_names)
  set.seed(20280000L + dataset * 10000L + method_seed_index * 100L + chain)
  initial_xi <- make_start(chain, dat$xi)
  neighbors <- make_neighbors(dat$X)
  fit <- DM_online_adaptive_VS_RJMCMC(
    X = dat$X, Y = dat$Y, s2 = rep(10, ncol(dat$Y)), r2 = rep(10, ncol(dat$Y)),
    a = 1, b = 9, g_initial = 1, c = 1, runs = runs, burn = burn, thin = thin,
    informed = spec$informed, adapt_g = spec$adapt, adaptation_start = 100L,
    adaptation_end = -1L, window_size = 25L, rm_scale = 0.5,
    decay_power = 0.6, min_g = 0.05, max_g = 10,
    neigh = if (spec$local_move) neighbors else NULL,
    use_move = spec$local_move, n_move_attempts = 1L, lambda_move = 1,
    adaptation_metric = "diversity",
    target_effective_neighbors = min(ncol(dat$X), max(5, 0.025 * ncol(dat$X))),
    adaptation_offset = 100, balance_directions = spec$balanced,
    uniform_mixture = spec$mixture, score_optimizer_maxit = score_maxit,
    coefficient_optimizer_maxit = coefficient_maxit, initial_xi = initial_xi,
    coefficient_warm_start = coefficient_warm_start,
    within_model_update = within_model_update,
    within_model_attempts = within_attempts, within_model_scale = within_scale,
    cache_score_fits = if (is.null(cache_override)) spec$cache else cache_override,
    discrete_proposal = spec$proposal,
    laplace_optimizer_maxit = spec$laplace_maxit
  )
  xi <- fit$xi
  model_size <- vapply(seq_len(dim(xi)[3]), function(s) sum(xi[, , s]), numeric(1))
  pips <- apply(xi, c(1, 2), mean)
  beta_draws <- full_beta_draws(fit)
  beta_mean <- apply(beta_draws, c(1, 2), mean)
  beta_lower <- apply(beta_draws, c(1, 2), quantile, probs = 0.025)
  beta_upper <- apply(beta_draws, c(1, 2), quantile, probs = 0.975)
  list(
    fit = fit, model_size = model_size, pips = pips,
    beta_mean = beta_mean, beta_lower = beta_lower, beta_upper = beta_upper,
    quality = quality
  )
}

record_fit <- function(result, dat, dataset, chain, method) {
  if (save_full_fits) {
    dm_validate_fit(result$fit, paste(dataset, method, chain))
    saveRDS(list(fit = result$fit, data = dat, dataset = dataset, chain = chain,
                 method = method, quality = result$quality, source_md5 = source_hashes,
                 settings = list(runs = runs, burn = burn, thin = thin,
                   score_optimizer_maxit = result$fit$score_optimizer_maxit,
                   coefficient_optimizer_maxit = coefficient_maxit)),
            file.path(out_dir, "full_fits", sprintf("dataset_%03d_%s_chain_%02d_%s.rds",
                       dataset, method, chain, result$quality)))
  }
  chain_id <<- chain_id + 1L
  fit <- result$fit
  model_size <- result$model_size
  chain_rows[[chain_id]] <<- data.frame(
    dataset = dataset, chain = chain, method = method, plugin_quality = result$quality,
    seconds = fit$total_seconds, iterations_per_second = fit$iterations_per_second,
    acceptance_flip = fit$acpt_rate_betaxi_global,
    acceptance_local = fit$acpt_rate_move_global,
    coefficient_warm_start = fit$coefficient_warm_start,
    within_model_update = fit$within_model_update,
    within_model_attempts = fit$within_model_attempts,
    within_model_scale = fit$within_model_scale,
    acceptance_within = if (sum(fit$within_model_proposed_by_category)>0)
      sum(fit$within_model_accepted_by_category)/sum(fit$within_model_proposed_by_category) else NA_real_,
    final_g_mean = mean(fit$g_final),
    mean_effective_neighbors = mean(fit$mean_proposal_effective_neighbors),
    mean_proposal_entropy = mean(fit$mean_proposal_entropy),
    mean_max_probability = mean(fit$mean_max_proposal_probability),
    ess_model_size = ess_basic(model_size),
    ess_model_size_per_second = ess_basic(model_size) / fit$total_seconds,
    bulk_ess_model_size = bulk_ess(model_size),
    tail_ess_model_size = tail_ess(model_size),
    ess_log_posterior = ess_basic(fit$log_posterior_trace),
    ess_log_posterior_per_second = ess_basic(fit$log_posterior_trace) / fit$total_seconds,
    mean_model_size = mean(model_size),
    pip_rmse = sqrt(mean((result$pips - dat$xi)^2)),
    mean_active_pip = mean(result$pips[dat$xi]),
    mean_inactive_pip = mean(result$pips[!dat$xi]),
    true_positive_rate_05 = mean(result$pips[dat$xi] >= 0.5),
    false_positive_rate_05 = mean(result$pips[!dat$xi] >= 0.5),
    active_coefficient_bias = mean(result$beta_mean[dat$xi] - dat$beta[dat$xi]),
    active_coefficient_rmse = sqrt(mean((result$beta_mean[dat$xi] - dat$beta[dat$xi])^2)),
    active_coefficient_coverage = mean(
      result$beta_lower[dat$xi] <= dat$beta[dat$xi] &
        result$beta_upper[dat$xi] >= dat$beta[dat$xi]
    ),
    score_optimizer_maxit = fit$score_optimizer_maxit,
    coefficient_optimizer_maxit = fit$coefficient_optimizer_maxit,
    cache_score_fits = fit$cache_score_fits,
    discrete_proposal = fit$discrete_proposal,
    laplace_optimizer_maxit = fit$laplace_optimizer_maxit,
    rhat_model_size = NA_real_, rhat_log_posterior = NA_real_
  )
  pip_id <<- pip_id + 1L
  pip_rows[[pip_id]] <<- data.frame(
    dataset = dataset, chain = chain, method = method, plugin_quality = result$quality,
    predictor = rep(seq_len(nrow(dat$xi)), times = ncol(dat$xi)),
    category = rep(seq_len(ncol(dat$xi)), each = nrow(dat$xi)),
    truth = as.vector(dat$xi), pip = as.vector(result$pips)
  )
  coefficient_id <<- coefficient_id + 1L
  coefficient_rows[[coefficient_id]] <<- data.frame(
    dataset = dataset, chain = chain, method = method, plugin_quality = result$quality,
    predictor = rep(seq_len(nrow(dat$beta)), times = ncol(dat$beta)),
    category = rep(seq_len(ncol(dat$beta)), each = nrow(dat$beta)),
    truth = as.vector(dat$beta), mean = as.vector(result$beta_mean),
    lower = as.vector(result$beta_lower), upper = as.vector(result$beta_upper),
    covered = as.vector(result$beta_lower <= dat$beta & result$beta_upper >= dat$beta)
  )
  keep <- unique(round(seq(1, length(model_size), length.out = min(500L, length(model_size)))))
  trace_id <<- trace_id + 1L
  trace_rows[[trace_id]] <<- data.frame(
    dataset = dataset, chain = chain, method = method, plugin_quality = result$quality,
    retained_iteration = keep, model_size = model_size[keep],
    log_posterior = fit$log_posterior_trace[keep]
  )
  adaptation_id <<- adaptation_id + 1L
  adaptation_rows[[adaptation_id]] <<- data.frame(
    dataset = dataset, chain = chain, method = method, plugin_quality = result$quality,
    retained_iteration = rep(keep, times = ncol(fit$g_trace)),
    category = rep(seq_len(ncol(fit$g_trace)), each = length(keep)),
    g = as.vector(fit$g_trace[keep, , drop = FALSE]),
    effective_neighbors = as.vector(fit$proposal_effective_neighbors_trace[keep, , drop = FALSE])
  )
}

if (stage %in% c("all", "core")) {
  for (dataset in dataset_ids) {
    dat <- simulate_dm(20281000L + dataset)
    first_chain_row <- chain_id + 1L
    size_traces <- setNames(vector("list", nrow(core_specs)), core_specs$method)
    size_traces <- lapply(size_traces, function(x) vector("list", chains))
    posterior_traces <- setNames(vector("list", nrow(core_specs)), core_specs$method)
    posterior_traces <- lapply(posterior_traces, function(x) vector("list", chains))
    for (spec_index in seq_len(nrow(core_specs))) {
      spec <- core_specs[spec_index, ]
      for (chain in seq_len(chains)) {
        result <- tryCatch(run_one(dat, dataset, chain, spec,
                                  score_maxit_core, cache_override = core_cache_override),
                           error = identity)
        if (inherits(result, "error")) {
          failure_id <- failure_id + 1L
          failures[[failure_id]] <- data.frame(
            dataset = dataset, chain = chain, method = spec$method,
            plugin_quality = "production", message = conditionMessage(result)
          )
          next
        }
        record_fit(result, dat, dataset, chain, spec$method)
        size_traces[[spec$method]][[chain]] <- result$model_size
        posterior_traces[[spec$method]][[chain]] <- result$fit$log_posterior_trace
        message(sprintf("DM core dataset=%d/%d method=%s chain=%d/%d",
                        dataset, datasets, spec$method, chain, chains))
      }
    }
    if (chain_id >= first_chain_row) {
      for (method in names(size_traces)) {
        available <- Filter(Negate(is.null), size_traces[[method]])
        available_posterior <- Filter(Negate(is.null), posterior_traces[[method]])
        rhat <- if (length(available) >= 2L) rank_split_rhat(available) else NA_real_
        posterior_rhat <- if (length(available_posterior) >= 2L) {
          rank_split_rhat(available_posterior)
        } else NA_real_
        indices <- which(vapply(chain_rows, function(x) {
          x$dataset == dataset && x$method == method && x$plugin_quality == "production"
        }, logical(1)))
        for (index in indices) {
          chain_rows[[index]]$rhat_model_size <- rhat
          chain_rows[[index]]$rhat_log_posterior <- posterior_rhat
        }
      }
    }
  }
}

if (stage %in% c("all", "plugin")) {
  dat <- simulate_dm(20289999L)
  spec <- core_specs[core_specs$method == "online_similarity_flip", ]
  quality_grid <- data.frame(
    quality = c("degraded", "default", "accurate"),
    score_maxit = c(1L, 1L, 10L), cache = c(FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  for (quality_index in seq_len(nrow(quality_grid))) {
    for (chain in seq_len(chains)) {
      quality <- quality_grid[quality_index, ]
      result <- tryCatch(
        run_one(dat, 0L, chain, spec, quality$score_maxit, quality$quality,
                quality$cache),
        error = identity
      )
      if (inherits(result, "error")) {
        failure_id <- failure_id + 1L
        failures[[failure_id]] <- data.frame(
          dataset = 0L, chain = chain, method = "online_similarity_flip",
          plugin_quality = quality$quality, message = conditionMessage(result)
        )
        next
      }
      record_fit(result, dat, 0L, chain, "online_similarity_flip")
      message(sprintf("DM plugin quality=%s chain=%d/%d", quality$quality, chain, chains))
    }
  }
}

bind_rows <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
chain_metrics <- bind_rows(chain_rows)
pips <- bind_rows(pip_rows)
coefficients <- bind_rows(coefficient_rows)
traces <- bind_rows(trace_rows)
adaptation <- bind_rows(adaptation_rows)
failure_table <- bind_rows(failures)

write.csv(chain_metrics, file.path(out_dir, "chain_metrics.csv"), row.names = FALSE)
write.csv(pips, file.path(out_dir, "pips.csv"), row.names = FALSE)
write.csv(coefficients, file.path(out_dir, "coefficient_summaries.csv"), row.names = FALSE)
write.csv(traces, file.path(out_dir, "model_size_traces_thinned.csv"), row.names = FALSE)
write.csv(adaptation, file.path(out_dir, "adaptive_trajectories_thinned.csv"), row.names = FALSE)
write.csv(failure_table, file.path(out_dir, "failures.csv"), row.names = FALSE)

if (nrow(chain_metrics)) {
  summary_columns <- setdiff(
    names(chain_metrics),
    c("dataset", "chain", "method", "plugin_quality", "cache_score_fits",
      "discrete_proposal")
  )
  method_summary <- aggregate(
    chain_metrics[summary_columns],
    by = list(method = chain_metrics$method,
              plugin_quality = chain_metrics$plugin_quality),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  write.csv(method_summary, file.path(out_dir, "overall_method_summary.csv"),
            row.names = FALSE)
}

if (nrow(pips)) {
  aggregated_pips <- aggregate(pip ~ dataset + method + plugin_quality + predictor + category + truth,
                               pips, mean)
  write.csv(aggregated_pips, file.path(out_dir, "aggregated_pips.csv"), row.names = FALSE)
  method_pairs <- merge(aggregated_pips, aggregated_pips,
                        by = c("dataset", "predictor", "category", "truth"),
                        suffixes = c("_a", "_b"))
  method_pairs <- method_pairs[method_pairs$method_a < method_pairs$method_b &
                                 method_pairs$plugin_quality_a == "production" &
                                 method_pairs$plugin_quality_b == "production", ]
  if (nrow(method_pairs)) {
    agreement <- aggregate((pip_a - pip_b)^2 ~ dataset + method_a + method_b,
                           method_pairs, function(x) sqrt(mean(x)))
    names(agreement)[4] <- "pip_rmse"
    write.csv(agreement, file.path(out_dir, "cross_kernel_pip_agreement.csv"), row.names = FALSE)
  }
}

writeLines(c(
  paste("mode:", mode), paste("stage:", stage),
  paste("dataset_ids:", paste(dataset_ids, collapse = ",")),
  paste("methods:", paste(core_specs$method, collapse = ",")),
  paste("chains:", chains), paste("runs:", runs), paste("burn:", burn),
  paste("thin:", thin), "n: 70", "J: 5", "P: 18", "predictor_correlation: AR(1), rho=0.65",
  paste("save_full_fits:", save_full_fits),
  paste("core_score_optimizer_maxit:", score_maxit_core),
  paste("core_cache_override:", cache_option),
  "current intercept covariance schedule: adapted during burn-in, then frozen",
  "production adaptation: category-wise online diminishing diversity adaptation",
  "adaptation gain: 0.5*(100+m)^(-0.6)",
  "effective-neighbor target: min(P,max(5,0.025*P))",
  "optimized proposal: balanced birth/death directions plus 2% direction-preserving uniform mixture",
  "continuous coefficient proposal: PMLE/Hessian fixed across uniform/similarity ablations",
  "Gagnon comparator: conditional posterior Laplace model weights with Barker balancing and matching posterior-mode/Hessian Gaussian RJ proposal",
  "local graph: four strongest empirical absolute correlations per predictor"
), file.path(out_dir, "manifest.txt"))

saveRDS(list(chain_metrics = chain_metrics, pips = pips, coefficients = coefficients,
             traces = traces, adaptation = adaptation, failures = failure_table),
        file.path(out_dir, "chain_results.rds"))

if (nrow(failure_table)) {
  stop(sprintf("%d chains failed; partial diagnostics were written to %s",
               nrow(failure_table), file.path(out_dir, "failures.csv")))
}

print(chain_metrics)
