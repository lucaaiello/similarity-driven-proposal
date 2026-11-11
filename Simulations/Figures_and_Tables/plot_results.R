experiment_folders <- c(exact="Exact_Gaussian",scaling="Gaussian_Scaling",swap="Correlated_Predictors",exact_score="Exact_Score_Sensitivity",correlated_score="Correlated_Score_Sensitivity")
# Plot saved results only: this script never launches a sampler.
# Rscript Simulations/Figures_and_Tables/plot_results.R
# Optional: --figure=gaussian|dm|all --gaussian-root=PATH --image-dir=PATH
# --application-root=PATH explicitly opts in to plotting saved application results.
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = "") {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 before running this script.")
suppressPackageStartupMessages(library(ggplot2))
selection <- option("figure", "all")
if (!selection %in% c("all", "gaussian", "dm")) stop("Invalid --figure")
output_root <- file.path(repo_root,"Simulations")
gaussian_root <- option("gaussian-root", file.path(repo_root,"Simulations"))
image_dir <- option("image-dir", file.path(repo_root,"reproduced/figures"))
application_root <- option("application-root")
dir.create(image_dir, recursive = TRUE, showWarnings = FALSE)

# Shared theme for main-text and supplementary panels. No plot titles or tags.
source(file.path(repo_root, "Simulations/Figures_and_Tables/plot_theme.R"))
method_labels <- c(uniform = "Uniform", f_fixed = "F, fixed", f_online = "F, adaptive",
                   bic_fixed = "BIC/LR, fixed", bic_online = "BIC/LR, adaptive",
                   locally_balanced = "Exact LB", bic_locally_balanced = "BIC/Laplace LB",
                   f_flip_only = "F: flip only", f_uniform_swap = "F: uniform swap",
                   f_informed_swap = "F: informed swap", bic_flip_only = "BIC/LR: flip only",
                   bic_uniform_swap = "BIC/LR: uniform swap", bic_informed_swap = "BIC/LR: informed swap",
                   uniform_flip = "Uniform", optimized_similarity_fixed = "Similarity, fixed",
                   online_similarity_flip = "Similarity, adaptive",
                   online_similarity_local = "Adaptive + swap", gagnon_laplace_rj = "Gagnon-Laplace")
score_order <- c("f", "bic", "marginal")
score_labels <- c(f = "F", bic = "BIC/LR", marginal = "Marginal")
exact_order <- c("uniform", "f_fixed", "f_online", "bic_fixed", "bic_online",
                 "locally_balanced", "bic_locally_balanced")
manifest <- list()
save_panel <- function(p, name, width = 4.8, height = 3.8) {
  # Padding is controlled by the discrete scale/ggplot, not manually shifted x limits.
  ggsave(file.path(image_dir, paste0(name, ".pdf")), p, device = grDevices::cairo_pdf,
         width = width, height = height, units = "in", bg = "white")
  ggsave(file.path(image_dir, paste0(name, ".png")), p,
         width = width, height = height, units = "in", dpi = 180, bg = "white")
  manifest[[length(manifest) + 1L]] <<- data.frame(file = paste0(name, ".pdf"), width, height)
}
read_result <- function(directory, filename) {
  path <- file.path(directory, filename)
  if (!file.exists(path)) stop("Missing result file: ", path)
  read.csv(path, check.names = FALSE)
}
aggregate_chains <- function(data, groups, variables) {
  stats::aggregate(data[variables], data[groups], mean, na.rm = TRUE)
}
summarize_grouped <- function(data, groups, variables) {
  pieces <- split(data, interaction(data[groups], drop = TRUE, lex.order = TRUE))
  do.call(rbind, lapply(pieces, function(z) {
    row <- z[1, groups, drop = FALSE]
    for (v in variables) {
      values <- z[[v]][is.finite(z[[v]])]
      row[[paste0(v, "_mean")]] <- if (length(values)) mean(values) else NA_real_
      row[[paste0(v, "_se")]] <- if (length(values) > 1) sd(values) / sqrt(length(values)) else NA_real_
    }
    row
  }))
}
metric_plot <- function(data, metric, order, ylabel, multiplier = 1, log_y = FALSE) {
  data$method <- factor(data$method, levels = order)
  data$value <- multiplier * data[[paste0(metric, "_mean")]]
  data$lower <- data$value - 1.96 * multiplier * data[[paste0(metric, "_se")]]
  data$upper <- data$value + 1.96 * multiplier * data[[paste0(metric, "_se")]]
  if (log_y && any(data$lower <= 0, na.rm = TRUE)) stop("Nonpositive interval on logarithmic axis")
  p <- ggplot(data, aes(method, value, colour = method)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.16, linewidth = 0.5) +
    geom_point(size = 2.4) + scale_colour_manual(values = palette, guide = "none") +
    scale_x_discrete(labels = method_labels, drop = FALSE) +
    labs(x = NULL, y = ylabel) + theme_manuscript(diagonal = TRUE)
  if (log_y) p <- p + scale_y_log10()
  p
}

if (selection != "dm") {
  # Intentionally default to the matched reruns, never silently to unmatched results.
  for (directory in c("exact", "scaling", "swap", "exact_score", "correlated_score")) {
    if (!file.exists(file.path(gaussian_root, experiment_folders[[directory]], "Results/provenance.rds")))
      stop("Missing matched Gaussian provenance: ", directory)
  }
  exact <- read_result(file.path(gaussian_root,"Exact_Gaussian/Results"), "chain_summary.csv")
  exact$robust_ess_second <- pmin(exact$ess_model_size_per_second, exact$ess_log_target_per_second)
  ed <- aggregate_chains(exact, c("n", "p", "dataset", "method"),
                         c("mean_abs_pip_error", "robust_ess_second"))
  es <- summarize_grouped(ed, "method", c("mean_abs_pip_error", "robust_ess_second"))
  save_panel(metric_plot(es, "mean_abs_pip_error", exact_order, "PIP MAE (percentage points)", 100),
             "simulation_exact_reference_accuracy_ggplot")
  save_panel(metric_plot(es, "robust_ess_second", exact_order, "Robust ESS/s", log_y = TRUE),
             "simulation_exact_reference_efficiency_ggplot")

  high <- read_result(file.path(gaussian_root,"Gaussian_Scaling/Results"), "chain_metrics.csv")
  high$robust_ess_second <- pmin(high$ess_model_size_per_second, high$ess_log_target_per_second)
  hd <- aggregate_chains(high, c("dimension", "dataset", "method"),
                         c("iterations_per_second", "robust_ess_second"))
  hs <- summarize_grouped(hd, c("dimension", "method"),
                          c("iterations_per_second", "robust_ess_second"))
  hs$method <- factor(hs$method, levels = exact_order)
  for (metric in c("iterations_per_second", "robust_ess_second")) {
    hs$value <- hs[[paste0(metric, "_mean")]]
    p <- ggplot(hs, aes(dimension, value, colour = method, group = method)) +
      geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
      scale_colour_manual(values = palette, labels = method_labels) +
      scale_x_log10(breaks = sort(unique(hs$dimension))) + scale_y_log10() +
      labs(x = "Number of predictors, P", y = if (metric == "iterations_per_second") "Transitions/s" else "Robust ESS/s") +
      theme_manuscript() + guides(colour = guide_legend(ncol = 2))
    save_panel(p, paste0("simulation_highdim_scaling_", if (metric == "iterations_per_second") "throughput" else "efficiency", "_ggplot"), height = 4.2)
  }

  local <- read_result(file.path(gaussian_root,"Correlated_Predictors/Results"), "chain_metrics.csv")
  local$robust_ess_second <- pmin(local$ess_model_size_per_second, local$ess_log_target / local$seconds)
  ld <- aggregate_chains(local, c("dataset", "method"), "robust_ess_second")
  blocks <- read_result(file.path(gaussian_root,"Correlated_Predictors/Results"), "block_metrics.csv")
  bd <- aggregate_chains(blocks, c("dataset", "method"), "switches_per_1000")
  ls <- summarize_grouped(merge(ld, bd, by = c("dataset", "method")), "method",
                          c("robust_ess_second", "switches_per_1000"))
  local_order <- c("f_flip_only", "f_uniform_swap", "f_informed_swap", "bic_flip_only", "bic_uniform_swap", "bic_informed_swap")
  save_panel(metric_plot(ls, "switches_per_1000", local_order, "Block-state changes / 1,000 iterations"),
             "simulation_correlated_local_moves_movement_ggplot")
  save_panel(metric_plot(ls, "robust_ess_second", local_order, "Robust ESS/s"),
             "simulation_correlated_local_moves_efficiency_ggplot")

  alignment <- read_result(file.path(gaussian_root,"Exact_Score_Sensitivity/Results"), "score_alignment.csv")
  alignment$score <- factor(alignment$score, levels = score_order)
  p <- ggplot(alignment, aes(score, posterior_weighted_rank_correlation, fill = score)) +
    geom_boxplot(width = 0.55, linewidth = 0.4, outlier.size = 1.3) +
    scale_fill_manual(values = palette, guide = "none") + scale_x_discrete(labels = score_labels) +
    labs(x = NULL, y = NULL) + theme_manuscript()
  save_panel(p, "simulation_auxiliary_score_alignment_ggplot", width = 3.6, height = 3.8)
  scores <- read_result(file.path(gaussian_root,"Exact_Score_Sensitivity/Results"), "chain_metrics.csv")
  scores$robust_ess_second <- pmin(scores$ess_model_size_per_second, scores$ess_log_target_per_second)
  sd <- aggregate_chains(scores, c("dataset", "score", "regime", "lambda_initial"), c("lambda_final", "robust_ess_second"))
  ss <- summarize_grouped(sd, c("score", "regime", "lambda_initial"), c("lambda_final", "robust_ess_second"))
  ss$score <- factor(ss$score, levels = score_order)
  fixed <- ss[ss$regime == "fixed", ]
  online <- ss[ss$regime == "online", ]
  p <- ggplot(fixed, aes(lambda_initial, robust_ess_second_mean, colour = score, group = score)) +
    geom_line(linewidth = 0.6) + geom_point(aes(shape = "Fixed"), size = 2) +
    geom_point(data = online, aes(x = lambda_final_mean, shape = "Adaptive"), size = 2.8) +
    scale_colour_manual(values = palette, labels = score_labels) +
    scale_shape_manual(values = c(Fixed = 1, Adaptive = 18)) +
    labs(x = expression(lambda), y = "Robust ESS/s") + theme_manuscript() +
    guides(colour = guide_legend(nrow = 1, order = 1, override.aes = list(shape = 1)),
           shape = guide_legend(nrow = 1, order = 2))
  save_panel(p, "simulation_auxiliary_score_tuning_ggplot", width = 3.6, height = 3.8)
  corr <- read_result(file.path(gaussian_root,"Correlated_Score_Sensitivity/Results"), "chain_metrics.csv")
  cd <- aggregate_chains(corr, c("dataset", "score"), c("mean_active_pip", "mean_surrogate_pip", "mean_noise_pip"))
  cs <- summarize_grouped(cd, "score", c("mean_active_pip", "mean_surrogate_pip", "mean_noise_pip"))
  pip <- do.call(rbind, lapply(c(Active = "mean_active_pip", Surrogate = "mean_surrogate_pip", Noise = "mean_noise_pip"),
                               function(metric) data.frame(score = cs$score, value = cs[[paste0(metric, "_mean")]], metric)))
  pip$group <- factor(pip$metric, levels = c("mean_active_pip", "mean_surrogate_pip", "mean_noise_pip"),
                      labels = c("Active", "Surrogate", "Noise"))
  pip$score <- factor(pip$score, levels = score_order)
  p <- ggplot(pip, aes(score, value, colour = group, shape = group)) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c(Active = "#0072B2", Surrogate = "#D55E00", Noise = "#666666")) +
    scale_shape_manual(values = c(Active = 16, Surrogate = 17, Noise = 15)) +
    scale_x_discrete(labels = score_labels, expand = expansion(add = 0.6)) +
    scale_y_continuous(limits = c(0, 1)) + labs(x = NULL, y = "Mean PIP") + theme_manuscript() +
    guides(colour = guide_legend(nrow = 1))
  save_panel(p, "simulation_auxiliary_score_pips_ggplot", width = 3.6, height = 3.8)
}

if (selection != "gaussian") {
  dm <- read_result(file.path(repo_root,"Simulations/Dirichlet_Multinomial/Comparison/Results/Combined"), "combined_chain_metrics.csv")
  dm$robust_ess_second <- pmin(dm$ess_model_size_per_second, dm$ess_log_posterior_per_second)
  dd <- aggregate_chains(dm, c("dataset", "method"), c("pip_rmse", "robust_ess_second"))
  ds <- summarize_grouped(dd, "method", c("pip_rmse", "robust_ess_second"))
  dm_order <- c("uniform_flip", "optimized_similarity_fixed", "online_similarity_flip", "online_similarity_local", "gagnon_laplace_rj")
  # Omit historical/internal implementations, as in the current manuscript.
  ds <- ds[ds$method %in% dm_order, ]
  save_panel(metric_plot(ds, "pip_rmse", dm_order, "PIP RMSE"), "simulation_dm_comparison_recovery_ggplot")
  save_panel(metric_plot(ds, "robust_ess_second", dm_order, "Robust ESS/s"), "simulation_dm_comparison_efficiency_ggplot")
}

# Optional application exports use exactly the same theme. No implicit use of pilots.
if (nzchar(application_root)) {
  pips <- read_result(application_root, "aggregated_pips.csv")
  traces <- read_result(application_root, "model_size_logposterior_traces.csv")
  for (method in intersect(c("online_similarity_flip", "online_similarity_local"), unique(pips$method))) {
    z <- pips[pips$method == method, ]
    z <- z[order(z$category, z$predictor), ]
    z$index <- seq_len(nrow(z))
    # Compute the Bayesian FDR rule from the supplied results, not an old fixed cutoff.
    ranked <- sort(z$mean_pip, decreasing = TRUE)
    eligible <- which(cumsum(1 - ranked) / seq_along(ranked) <= 0.05)
    cutoff <- if (length(eligible)) ranked[max(eligible)] else Inf
    z$selected <- z$mean_pip >= cutoff
    suffix <- if (method == "online_similarity_local") "_move" else ""
    p <- ggplot(z, aes(index, mean_pip)) + geom_linerange(aes(ymin = 0, ymax = mean_pip), linewidth = 0.25, colour = "#666666") +
      geom_point(data = z[z$selected, ], colour = "#D55E00", size = 1.4) +
      labs(x = "Predictor-category index", y = "PIP") + theme_manuscript()
    save_panel(p, paste0("PIPs_adaptive", suffix, "_ggplot"), width = 7.2, height = 3.3)
    tr <- traces[traces$method == method, ]
    tr$chain <- factor(tr$chain)
    p <- ggplot(tr, aes(retained_iteration, model_size, group = chain)) + geom_line(linewidth = 0.35, colour = "#0072B2") +
      facet_wrap(~chain, ncol = 1, labeller = label_both) + labs(x = "Retained draw index", y = "Active associations") + theme_manuscript()
    save_panel(p, paste0("traceplot_adaptive", suffix, "_ggplot"), height = 4.5)
    ac <- do.call(rbind, lapply(split(tr, tr$chain), function(z) {
      z <- z[order(z$retained_iteration), ]
      a <- stats::acf(z$model_size, plot = FALSE, lag.max = min(100, nrow(z) - 1))
      data.frame(chain = z$chain[1], lag = as.vector(a$lag), value = as.vector(a$acf))
    }))
    p <- ggplot(ac, aes(lag, value)) + geom_segment(aes(xend = lag, yend = 0), linewidth = 0.35, colour = "#0072B2") +
      facet_wrap(~chain, ncol = 1, labeller = label_both) + labs(x = "Lag (retained draws)", y = "Autocorrelation") + theme_manuscript()
    save_panel(p, paste0("autocorrelation_adaptive", suffix, "_ggplot"), height = 4.5)
  }
}
write.csv(do.call(rbind, manifest), file.path(image_dir, "ggplot_panel_manifest.csv"), row.names = FALSE)
message("Saved separate ggplot panels to: ", image_dir)
