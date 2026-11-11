# Read-only by default. --output=<fresh directory> writes diagnostic summaries.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = "") {
  x <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(x)) sub(paste0("^--", name, "="), "", x[1]) else default
}
inputs <- strsplit(option("inputs", paste(
  "New_simulations/outputs/dm_online_adaptive_primary_20datasets/dm_online_adaptive_ablation_results.rds",
  "New_simulations/outputs/dm_gagnon_laplace_20datasets/dm_online_adaptive_ablation_results.rds",
  sep = ";")), ";", fixed = TRUE)[[1]]
out <- option("output")
pip_threshold <- as.numeric(option("pip-range-threshold", "0.10"))
if (!is.finite(pip_threshold) || pip_threshold <= 0) stop("Invalid PIP-range threshold")
objects <- lapply(inputs, readRDS)
methods <- c("uniform_flip", "optimized_similarity_fixed", "online_similarity_flip",
             "online_similarity_local", "gagnon_laplace_rj")
labels <- c("Uniform", "Similarity, fixed", "Similarity, adaptive",
            "Adaptive plus swap", "Gagnon--Laplace")
metric_columns <- c("dataset", "chain", "method", "plugin_quality", "rhat_model_size", "rhat_log_posterior")
metrics <- do.call(rbind, lapply(objects, function(x) x$chain_metrics[metric_columns]))
pips <- do.call(rbind, lapply(objects, `[[`, "pips"))
metrics <- metrics[metrics$method %in% methods & metrics$plugin_quality == "production", ]
pips <- pips[pips$method %in% methods & pips$plugin_quality == "production", ]
if (anyDuplicated(metrics[c("dataset", "method", "chain")])) stop("Duplicate chain records across inputs")
if (anyDuplicated(pips[c("dataset", "method", "chain", "predictor", "category")])) stop("Duplicate PIP records")
keys <- unique(metrics[c("dataset", "method")])
cases <- lapply(seq_len(nrow(keys)), function(i) {
  k <- keys[i, ]; z <- metrics[metrics$dataset == k$dataset & metrics$method == k$method, ]
  p <- pips[pips$dataset == k$dataset & pips$method == k$method, ]
  one <- function(name) {
    x <- unique(z[[name]])
    if (length(x) != 1L) stop("Inconsistent stored dataset-level Rhat")
    x
  }
  ranges <- if (nrow(p)) aggregate(pip ~ predictor + category, p,
                                 function(x) diff(range(x)))$pip else NA_real_
  rs <- one("rhat_model_size"); rl <- one("rhat_log_posterior")
  missing <- !is.finite(rs) || !is.finite(rl) || length(unique(z$chain)) < 3L ||
             !length(ranges) || any(!is.finite(ranges))
  mean_range <- mean(ranges); max_range <- max(ranges)
  data.frame(dataset = k$dataset, method = k$method, chains = nrow(z),
    rhat_model_size = rs, rhat_log_posterior = rl,
    mean_pip_range = mean_range, max_pip_range = max_range,
    missing_diagnostics = missing,
    above_1_01 = any(c(rs, rl) > 1.01, na.rm = TRUE),
    above_1_05 = any(c(rs, rl) > 1.05, na.rm = TRUE),
    pip_range_flag = is.finite(max_range) && max_range > pip_threshold)
})
cases <- do.call(rbind, cases)
cases$screen_flag <- with(cases, missing_diagnostics | above_1_01 | pip_range_flag)
triple <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, 3))
  c(median(x), unname(quantile(x, .75)), max(x))
}
summary <- do.call(rbind, lapply(seq_along(methods), function(i) {
  z <- cases[cases$method == methods[i], ]
  s <- triple(z$rhat_model_size); l <- triple(z$rhat_log_posterior)
  data.frame(method = methods[i], label = labels[i], datasets = nrow(z),
    size_median = s[1], size_q75 = s[2], size_max = s[3],
    log_median = l[1], log_q75 = l[2], log_max = l[3],
    n_above_1_01 = sum(z$above_1_01), n_above_1_05 = sum(z$above_1_05),
    n_missing = sum(z$missing_diagnostics),
    pip_range_median = median(z$mean_pip_range), pip_range_max = max(z$mean_pip_range),
    n_pip_flags = sum(z$pip_range_flag))
}))
print(summary, row.names = FALSE)
cat("\nScreened cases (screening is not proof of nonconvergence):\n")
print(cases[cases$screen_flag, ], row.names = FALSE)
fmt <- function(x, digits = 3L) formatC(x, format = "f", digits = digits)
latex <- vapply(seq_len(nrow(summary)), function(i) {
  z <- summary[i, ]
  paste(z$label,
        paste(fmt(c(z$size_median, z$size_q75, z$size_max)), collapse = "/"),
        paste(fmt(c(z$log_median, z$log_q75, z$log_max)), collapse = "/"),
        z$n_above_1_01, z$n_above_1_05, z$n_missing,
        paste(fmt(c(z$pip_range_median, z$pip_range_max), 4L), collapse = "/"),
        sep = " & ")
}, character(1))
latex <- paste0(latex, " \\\\")
cat("\nLaTeX rows:\n", paste(latex, collapse = "\n"), "\n")
if (nzchar(out)) {
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a fresh output directory")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(summary, file.path(out, "diagnostic_summary.csv"), row.names = FALSE)
  write.csv(cases, file.path(out, "dataset_diagnostics.csv"), row.names = FALSE)
  write.csv(cases[cases$screen_flag, ], file.path(out, "problematic_cases.csv"), row.names = FALSE)
  writeLines(latex, file.path(out, "table_rows.tex"))
  hashes <- tools::md5sum(inputs)
  write.csv(data.frame(file = names(hashes), md5 = unname(hashes)),
            file.path(out, "input_md5.csv"), row.names = FALSE)
}
