args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "New_simulations/make_manuscript_simulation_figures.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
dm_only <- "--figure=dm" %in% commandArgs(trailingOnly = TRUE)
gaussian_only <- "--figure=gaussian" %in% commandArgs(trailingOnly = TRUE)
figure_args <- commandArgs(trailingOnly=TRUE)
figure_arg <- function(name,default) {
  z <- grep(paste0("^--",name,"="),figure_args,value=TRUE)
  if(length(z))sub(paste0("^--",name,"="),"",z[[1]])else default
}
gaussian_root <- figure_arg("gaussian-root","")
# The manuscript uses the matched reruns. Earlier development output trees
# are archived and are no longer fallback inputs to this plotting script.
gaussian_read_root <- if(nzchar(gaussian_root)) gaussian_root else
  file.path(repo_root,"New_simulations","outputs","matched_gaussian_20261007_163020")
corrected_directories <- c(
  production_all_methods_exact_grid_rank_one_20datasets="exact",
  highdim_all_methods_full="scaling",correlated_surrogate_local_move_full="swap",
  auxiliary_score_exact_sensitivity_full="exact_score",
  auxiliary_score_correlated_full="correlated_score")

output_root <- file.path(repo_root, "New_simulations", "outputs")
image_dir <- file.path(repo_root, "JCGS_submission", "Images")
summary_dir <- file.path(output_root, "manuscript_simulation_summaries")
if(nzchar(gaussian_root)) {
  image_dir <- file.path(gaussian_root,"Images")
  summary_dir <- file.path(gaussian_root,"manuscript_simulation_summaries")
  if(!gaussian_only)stop("Use --figure=gaussian with --gaussian-root to avoid mixing new Gaussian and old DM results.")
  for(directory in unname(corrected_directories))
    if(!file.exists(file.path(gaussian_root,directory,"revision_provenance.rds")))
      stop("Missing corrected experiment: ",directory)
}
image_dir <- figure_arg("image-dir", image_dir)
summary_dir <- figure_arg("summary-dir", summary_dir)
dir.create(image_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

read_output <- function(directory, file) {
  if(directory%in%names(corrected_directories))
    return(read.csv(file.path(gaussian_read_root,corrected_directories[[directory]],file),check.names=FALSE))
  read.csv(file.path(output_root, directory, file), check.names = FALSE)
}

se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))

summarize_grouped <- function(data, groups, variables) {
  keys <- interaction(data[groups], drop = TRUE, lex.order = TRUE)
  pieces <- split(data, keys)
  rows <- lapply(pieces, function(z) {
    identifiers <- z[1, groups, drop = FALSE]
    values <- unlist(lapply(variables, function(v) {
      x <- z[[v]]
      c(mean = mean(x, na.rm = TRUE), se = se(x))
    }))
    names(values) <- unlist(lapply(variables, function(v) paste0(v, c("_mean", "_se"))))
    cbind(identifiers, as.data.frame(as.list(values), check.names = FALSE))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

aggregate_chains <- function(data, groups, variables) {
  stats::aggregate(data[variables], data[groups], mean, na.rm = TRUE)
}

method_labels <- c(
  uniform = "Uniform",
  f_fixed = "F, fixed",
  f_online = "F, adaptive",
  bic_fixed = "BIC/LR, fixed",
  bic_online = "BIC/LR, adaptive",
  locally_balanced = "Exact LB",
  bic_locally_balanced = "BIC/Laplace LB",
  f_flip_only = "F: flip only",
  f_uniform_swap = "F: uniform swap",
  f_informed_swap = "F: informed swap",
  bic_flip_only = "BIC/LR: flip only",
  bic_uniform_swap = "BIC/LR: uniform swap",
  bic_informed_swap = "BIC/LR: informed swap",
  uniform_flip = "Uniform",
  legacy_similarity_flip = "Legacy similarity",
  optimized_similarity_fixed = "Similarity, fixed",
  online_similarity_flip = "Similarity, adaptive",
  online_similarity_local = "Adaptive + swap",
  gagnon_laplace_rj = "Gagnon--Laplace"
)

palette <- c(
  uniform = "#666666", f_fixed = "#56B4E9", f_online = "#0072B2",
  bic_fixed = "#E69F00", bic_online = "#D55E00",
  locally_balanced = "#009E73", bic_locally_balanced = "#CC79A7",
  f_flip_only = "#56B4E9", f_uniform_swap = "#0072B2", f_informed_swap = "#003F6B",
  bic_flip_only = "#E69F00", bic_uniform_swap = "#D55E00", bic_informed_swap = "#8C2D04",
  uniform_flip = "#666666", legacy_similarity_flip = "#999999",
  optimized_similarity_fixed = "#E69F00", online_similarity_flip = "#0072B2",
  online_similarity_local = "#009E73", gagnon_laplace_rj = "#CC79A7",
  f = "#0072B2", bic = "#D55E00", marginal = "#009E73"
)

draw_ci <- function(x, y, yse, col, pch = 19, horizontal = FALSE) {
  if (horizontal) {
    segments(y - 1.96 * yse, x, y + 1.96 * yse, x, col = col, lwd = 1.5)
    points(y, x, pch = pch, col = col, bg = col, cex = 1.05)
  } else {
    segments(x, y - 1.96 * yse, x, y + 1.96 * yse, col = col, lwd = 1.5)
    points(x, y, pch = pch, col = col, bg = col, cex = 1.05)
  }
}

if (!dm_only) {
## Figure 1 and exact-reference table ---------------------------------------
exact <- read_output("production_all_methods_exact_grid_rank_one_20datasets", "chain_summary.csv")
exact$robust_ess_second <- pmin(exact$ess_model_size_per_second, exact$ess_log_target_per_second)
exact_dataset <- aggregate_chains(
  exact,
  c("n", "p", "dataset", "method"),
  c("mean_abs_pip_error", "total_variation_distance", "robust_ess_second", "acceptance_rate",
    "final_lambda", "mean_last_quarter_effective_neighbors")
)
exact_summary <- summarize_grouped(
  exact_dataset, "method",
  c("mean_abs_pip_error", "total_variation_distance", "robust_ess_second", "acceptance_rate",
    "final_lambda", "mean_last_quarter_effective_neighbors")
)
write.csv(exact_summary, file.path(summary_dir, "table_exact_reference.csv"), row.names = FALSE)

exact_order <- c("uniform", "f_fixed", "f_online", "bic_fixed", "bic_online",
                 "locally_balanced", "bic_locally_balanced")
es <- exact_summary[match(exact_order, exact_summary$method), ]

pdf(file.path(image_dir, "simulation_exact_reference.pdf"), width = 9.2, height = 4.8,
    family = "sans", pointsize = 10)
par(mfrow = c(1, 2), mar = c(7.1, 4.2, 2.0, 0.8), oma = c(0, 0, 0, 0))
x <- seq_along(exact_order)
plot(x, 100 * es$mean_abs_pip_error_mean, type = "n", xaxt = "n",
     xlab = "", ylab = "Mean absolute PIP error (percentage points)",
     ylim = range(100 * c(es$mean_abs_pip_error_mean - 1.96 * es$mean_abs_pip_error_se,
                         es$mean_abs_pip_error_mean + 1.96 * es$mean_abs_pip_error_se)),
     main = "(a) Exact-posterior accuracy")
axis(1, at = x, labels = method_labels[exact_order], las = 2, cex.axis = 0.82)
draw_ci(x, 100 * es$mean_abs_pip_error_mean, 100 * es$mean_abs_pip_error_se,
        palette[exact_order])
grid(nx = NA, ny = NULL, col = "#E5E5E5")

plot(x, es$robust_ess_second_mean, type = "n", xaxt = "n", log = "y",
     xlab = "", ylab = "Robust ESS per second (log scale)",
     ylim = range(c(es$robust_ess_second_mean - 1.96 * es$robust_ess_second_se,
                    es$robust_ess_second_mean + 1.96 * es$robust_ess_second_se)),
     main = "(b) Sampling efficiency")
axis(1, at = x, labels = method_labels[exact_order], las = 2, cex.axis = 0.82)
draw_ci(x, es$robust_ess_second_mean, es$robust_ess_second_se, palette[exact_order])
grid(nx = NA, ny = NULL, col = "#E5E5E5")
dev.off()

## Figure 2 and high-dimensional table -------------------------------------
high <- read_output("highdim_all_methods_full", "chain_metrics.csv")
high$robust_ess_second <- pmin(high$ess_model_size_per_second, high$ess_log_target_per_second)
high_dataset <- aggregate_chains(
  high, c("dimension", "dataset", "method"),
  c("iterations_per_second", "robust_ess_second", "top_signal_recall", "mean_inactive_pip")
)
high_summary <- summarize_grouped(
  high_dataset, c("dimension", "method"),
  c("iterations_per_second", "robust_ess_second", "top_signal_recall", "mean_inactive_pip")
)
exponents <- read_output("highdim_all_methods_full", "scaling_exponents.csv")
high_table <- merge(high_summary[high_summary$dimension == 1000, ], exponents, by = "method", all.x = TRUE)
write.csv(high_table, file.path(summary_dir, "table_highdim_p1000.csv"), row.names = FALSE)

pdf(file.path(image_dir, "simulation_highdim_scaling.pdf"), width = 9.2, height = 4.5,
    family = "sans", pointsize = 10)
par(mfrow = c(1, 2), mar = c(4.4, 4.4, 2.0, 1.0))
dims <- sort(unique(high_summary$dimension))
for (metric in c("iterations_per_second_mean", "robust_ess_second_mean")) {
  yr <- range(high_summary[[metric]], finite = TRUE)
  ttl <- if (metric == "iterations_per_second_mean") "(a) Transition throughput" else "(b) Sampling efficiency"
  ylab <- if (metric == "iterations_per_second_mean") "Iterations per second" else "Robust ESS per second"
  plot(dims, dims, type = "n", log = "xy", ylim = yr, xlab = "Number of predictors, P",
       ylab = ylab, main = ttl)
  for (m in exact_order) {
    z <- high_summary[high_summary$method == m, ]
    z <- z[order(z$dimension), ]
    lines(z$dimension, z[[metric]], type = "b", pch = 19, lwd = 1.35,
          col = palette[m], cex = 0.75)
  }
  grid(col = "#E5E5E5")
}
legend("bottomleft", legend = method_labels[exact_order], col = palette[exact_order],
       lty = 1, pch = 19, cex = 0.66, bty = "n", ncol = 2)
dev.off()

## Figure 3 and local-move table -------------------------------------------
local <- read_output("correlated_surrogate_local_move_full", "chain_metrics.csv")
local$robust_ess_second <- pmin(local$ess_model_size_per_second,
                                local$ess_log_target / local$seconds)
local_dataset <- aggregate_chains(
  local, c("dataset", "method"),
  c("robust_ess_second", "accepted_swaps_per_iteration", "multi_component_jump_rate",
    "rhat_model_size")
)
blocks <- read_output("correlated_surrogate_local_move_full", "block_metrics.csv")
block_dataset <- aggregate_chains(blocks, c("dataset", "method"),
                                  c("switches_per_1000", "true_pip", "surrogate_pip_sum"))
local_dataset <- merge(local_dataset, block_dataset, by = c("dataset", "method"), all.x = TRUE)
local_summary <- summarize_grouped(
  local_dataset, "method",
  c("switches_per_1000", "robust_ess_second", "multi_component_jump_rate",
    "true_pip", "surrogate_pip_sum", "rhat_model_size")
)
write.csv(local_summary, file.path(summary_dir, "table_correlated_local_moves.csv"), row.names = FALSE)

local_order <- c("f_flip_only", "f_uniform_swap", "f_informed_swap",
                 "bic_flip_only", "bic_uniform_swap", "bic_informed_swap")
lsum <- local_summary[match(local_order, local_summary$method), ]

pdf(file.path(image_dir, "simulation_correlated_local_moves.pdf"), width = 9.2, height = 4.8,
    family = "sans", pointsize = 10)
par(mfrow = c(1, 2), mar = c(7.3, 4.3, 2.0, 0.8))
x <- seq_along(local_order)
plot(x, lsum$switches_per_1000_mean, type = "n", xaxt = "n", xlab = "",
     ylab = "Within-block switches per 1,000 iterations",
     ylim = range(c(lsum$switches_per_1000_mean - 1.96 * lsum$switches_per_1000_se,
                    lsum$switches_per_1000_mean + 1.96 * lsum$switches_per_1000_se)),
     main = "(a) Correlated-substitute exploration")
axis(1, at = x, labels = method_labels[local_order], las = 2, cex.axis = 0.8)
draw_ci(x, lsum$switches_per_1000_mean, lsum$switches_per_1000_se, palette[local_order])
grid(nx = NA, ny = NULL, col = "#E5E5E5")

plot(x, lsum$robust_ess_second_mean, type = "n", xaxt = "n", xlab = "",
     ylab = "Robust ESS per second",
     ylim = range(c(lsum$robust_ess_second_mean - 1.96 * lsum$robust_ess_second_se,
                    lsum$robust_ess_second_mean + 1.96 * lsum$robust_ess_second_se)),
     main = "(b) Overall sampling efficiency")
axis(1, at = x, labels = method_labels[local_order], las = 2, cex.axis = 0.8)
draw_ci(x, lsum$robust_ess_second_mean, lsum$robust_ess_second_se, palette[local_order])
grid(nx = NA, ny = NULL, col = "#E5E5E5")
dev.off()

## Figure 4 and auxiliary-score table --------------------------------------
alignment <- read_output("auxiliary_score_exact_sensitivity_full", "score_alignment.csv")
alignment_summary <- summarize_grouped(
  alignment, "score", c("posterior_weighted_rank_correlation", "unweighted_rank_correlation")
)
score_exact <- read_output("auxiliary_score_exact_sensitivity_full", "chain_metrics.csv")
score_exact$robust_ess_second <- pmin(score_exact$ess_model_size_per_second,
                                      score_exact$ess_log_target_per_second)
score_exact_dataset <- aggregate_chains(
  score_exact, c("dataset", "score", "regime", "lambda_initial"),
  c("lambda_final", "robust_ess_second", "mean_pip_error", "mean_effective_neighbors")
)
score_exact_summary <- summarize_grouped(
  score_exact_dataset, c("score", "regime", "lambda_initial"),
  c("lambda_final", "robust_ess_second", "mean_pip_error", "mean_effective_neighbors")
)
score_corr <- read_output("auxiliary_score_correlated_full", "chain_metrics.csv")
score_corr$robust_ess_second <- pmin(score_corr$ess_model_size_per_second,
                                     score_corr$ess_log_target_per_second)
score_corr_dataset <- aggregate_chains(
  score_corr, c("dataset", "score"),
  c("robust_ess_second", "mean_active_pip", "mean_surrogate_pip", "mean_noise_pip",
    "lambda_final", "mean_effective_neighbors")
)
score_corr_summary <- summarize_grouped(
  score_corr_dataset, "score",
  c("robust_ess_second", "mean_active_pip", "mean_surrogate_pip", "mean_noise_pip",
    "lambda_final", "mean_effective_neighbors")
)
score_table <- merge(alignment_summary, score_corr_summary, by = "score", all = TRUE)
write.csv(score_table, file.path(summary_dir, "table_auxiliary_score.csv"), row.names = FALSE)

score_order <- c("f", "bic", "marginal")
pdf(file.path(image_dir, "simulation_auxiliary_score.pdf"), width = 10.2, height = 4.5,
    family = "sans", pointsize = 10)
par(mfrow = c(1, 3), mar = c(4.8, 4.2, 2.0, 0.7))

abox <- split(alignment$posterior_weighted_rank_correlation,
              factor(alignment$score, levels = score_order))
boxplot(abox, names = c("F", "BIC/LR", "Marginal"), col = palette[score_order], border = "#444444",
        ylab = "Posterior-weighted rank correlation", xlab = "Auxiliary score",
        main = "(a) Score--target alignment")
grid(nx = NA, ny = NULL, col = "#E5E5E5")

fixed <- score_exact_summary[score_exact_summary$regime == "fixed", ]
yr <- range(fixed$robust_ess_second_mean, finite = TRUE)
plot(range(fixed$lambda_initial), yr, type = "n", xlab = expression(lambda),
     ylab = "Robust ESS per second", main = "(b) Fixed and adaptive tuning")
for (s in score_order) {
  z <- fixed[fixed$score == s, ]
  z <- z[order(z$lambda_initial), ]
  lines(z$lambda_initial, z$robust_ess_second_mean, type = "b", pch = 1,
        col = palette[s], lwd = 1.35)
  online <- score_exact_summary[score_exact_summary$score == s & score_exact_summary$regime == "online", ]
  points(online$lambda_final_mean, online$robust_ess_second_mean, pch = 18,
         col = palette[s], cex = 1.25)
}
grid(col = "#E5E5E5")
legend("topright", legend = c("F", "BIC/LR", "Marginal", "Adaptive result"),
       col = c(palette[score_order], "#222222"), pch = c(1, 1, 1, 18), lty = c(1, 1, 1, NA),
       bty = "n", cex = 0.72)

cs <- score_corr_summary[match(score_order, score_corr_summary$score), ]
x <- seq_along(score_order)
plot(x, cs$mean_active_pip_mean, type = "n", xaxt = "n", ylim = c(0, 1), xlab = "",
     ylab = "Mean posterior inclusion probability", main = "(c) Correlated design")
axis(1, at = x, labels = c("F", "BIC/LR", "Marginal"))
points(x - 0.12, cs$mean_active_pip_mean, pch = 19, col = "#0072B2")
points(x, cs$mean_surrogate_pip_mean, pch = 17, col = "#D55E00")
points(x + 0.12, cs$mean_noise_pip_mean, pch = 15, col = "#666666")
grid(nx = NA, ny = NULL, col = "#E5E5E5")
legend("topright", legend = c("Active", "Correlated surrogate", "Noise"),
       col = c("#0072B2", "#D55E00", "#666666"), pch = c(19, 17, 15), bty = "n", cex = 0.72)
dev.off()

}

## Figure 5 and Dirichlet--multinomial table --------------------------------
if (!gaussian_only) {
dm <- read_output("dm_gagnon_comparison", "combined_chain_metrics.csv")
dm$robust_ess_second <- pmin(dm$ess_model_size_per_second, dm$ess_log_posterior_per_second)
dm_dataset <- aggregate_chains(
  dm, c("dataset", "method"),
  c("pip_rmse", "robust_ess_second", "true_positive_rate_05", "false_positive_rate_05",
    "active_coefficient_rmse", "active_coefficient_coverage", "seconds",
    "rhat_model_size", "rhat_log_posterior")
)
dm_summary <- summarize_grouped(
  dm_dataset, "method",
  c("pip_rmse", "robust_ess_second", "true_positive_rate_05", "false_positive_rate_05",
    "active_coefficient_rmse", "active_coefficient_coverage", "seconds",
    "rhat_model_size", "rhat_log_posterior")
)
write.csv(dm_summary, file.path(summary_dir, "table_dm_comparison.csv"), row.names = FALSE)

dm_order <- c("uniform_flip", "optimized_similarity_fixed",
              "online_similarity_flip", "online_similarity_local", "gagnon_laplace_rj")
ds <- dm_summary[match(dm_order, dm_summary$method), ]

pdf(file.path(image_dir, "simulation_dm_comparison.pdf"), width = 9.2, height = 4.8,
    family = "sans", pointsize = 10)
par(mfrow = c(1, 2), mar = c(7.2, 4.3, 2.0, 0.8))
x <- seq_along(dm_order)
plot(x, ds$pip_rmse_mean, type = "n", xaxt = "n", xlab = "", ylab = "PIP RMSE",
     ylim = range(c(ds$pip_rmse_mean - 1.96 * ds$pip_rmse_se,
                    ds$pip_rmse_mean + 1.96 * ds$pip_rmse_se)),
     main = "(a) Variable-selection accuracy")
axis(1, at = x, labels = method_labels[dm_order], las = 2, cex.axis = 0.8)
draw_ci(x, ds$pip_rmse_mean, ds$pip_rmse_se, palette[dm_order])
grid(nx = NA, ny = NULL, col = "#E5E5E5")

plot(x, ds$robust_ess_second_mean, type = "n", xaxt = "n", xlab = "",
     ylab = "Robust ESS per second",
     ylim = range(c(ds$robust_ess_second_mean - 1.96 * ds$robust_ess_second_se,
                    ds$robust_ess_second_mean + 1.96 * ds$robust_ess_second_se)),
     main = "(b) Sampling efficiency")
axis(1, at = x, labels = method_labels[dm_order], las = 2, cex.axis = 0.8)
draw_ci(x, ds$robust_ess_second_mean, ds$robust_ess_second_se, palette[dm_order])
grid(nx = NA, ny = NULL, col = "#E5E5E5")
dev.off()

}
message("Wrote manuscript figures to: ", image_dir)
message("Wrote derived table summaries to: ", summary_dir)
