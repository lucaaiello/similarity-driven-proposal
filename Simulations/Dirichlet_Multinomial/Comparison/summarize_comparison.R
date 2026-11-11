args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

baseline_dir <- arg_value(
  "baseline",
  file.path("Simulations/Dirichlet_Multinomial/Comparison/Results/Similarity_Methods")
)
gagnon_dir <- arg_value(
  "gagnon",
  file.path("Simulations/Dirichlet_Multinomial/Comparison/Results/Gagnon_Laplace")
)
out_dir <- arg_value(
  "output",
  file.path("Simulations/Dirichlet_Multinomial/Comparison/Results/Combined")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(directory, filename) {
  path <- file.path(directory, filename)
  if (!file.exists(path)) stop("missing required file: ", path)
  read.csv(path, stringsAsFactors = FALSE)
}

baseline_chains <- read_required(baseline_dir, "chain_metrics.csv")
gagnon_chains <- read_required(gagnon_dir, "chain_metrics.csv")
baseline_pips <- read_required(baseline_dir, "pips.csv")
gagnon_pips <- read_required(gagnon_dir, "pips.csv")

if (!identical(sort(unique(baseline_chains$dataset)),
               sort(unique(gagnon_chains$dataset)))) {
  stop("baseline and Gagnon runs do not contain the same dataset IDs")
}
if (!all(gagnon_chains$method == "gagnon_laplace_rj")) {
  stop("the Gagnon directory contains unexpected methods")
}

common_columns <- intersect(names(baseline_chains), names(gagnon_chains))
chains <- rbind(
  baseline_chains[, common_columns, drop = FALSE],
  gagnon_chains[, common_columns, drop = FALSE]
)
numeric_metrics <- names(chains)[vapply(chains, is.numeric, logical(1))]
numeric_metrics <- setdiff(numeric_metrics, c("dataset", "chain"))
method_summary <- aggregate(
  chains[numeric_metrics], by = list(method = chains$method),
  FUN = function(x) mean(x, na.rm = TRUE)
)

common_pip_columns <- intersect(names(baseline_pips), names(gagnon_pips))
pips <- rbind(
  baseline_pips[, common_pip_columns, drop = FALSE],
  gagnon_pips[, common_pip_columns, drop = FALSE]
)
aggregated_pips <- aggregate(
  pip ~ dataset + method + predictor + category + truth, pips, mean
)
gagnon_pip <- aggregated_pips[aggregated_pips$method == "gagnon_laplace_rj", ]
comparators <- aggregated_pips[aggregated_pips$method != "gagnon_laplace_rj", ]
agreement <- merge(
  comparators, gagnon_pip,
  by = c("dataset", "predictor", "category", "truth"),
  suffixes = c("_comparator", "_gagnon")
)
agreement_summary <- aggregate(
  (pip_comparator - pip_gagnon)^2 ~ dataset + method_comparator,
  agreement, function(x) sqrt(mean(x))
)
names(agreement_summary)[3] <- "pip_rmse_vs_gagnon"

gagnon_timing <- gagnon_chains[, c("dataset", "chain", "seconds",
                                    "ess_model_size_per_second",
                                    "ess_log_posterior_per_second")]
names(gagnon_timing)[3:5] <- paste0(names(gagnon_timing)[3:5], "_gagnon")
baseline_timing <- baseline_chains[, c("dataset", "chain", "method", "seconds",
                                       "ess_model_size_per_second",
                                       "ess_log_posterior_per_second")]
timing <- merge(baseline_timing, gagnon_timing, by = c("dataset", "chain"))
timing$runtime_ratio_gagnon_to_comparator <- timing$seconds_gagnon / timing$seconds
timing$ess_model_size_per_second_ratio_gagnon_to_comparator <-
  timing$ess_model_size_per_second_gagnon / timing$ess_model_size_per_second
timing$ess_log_posterior_per_second_ratio_gagnon_to_comparator <-
  timing$ess_log_posterior_per_second_gagnon / timing$ess_log_posterior_per_second

write.csv(chains, file.path(out_dir, "combined_chain_metrics.csv"), row.names = FALSE)
write.csv(method_summary, file.path(out_dir, "combined_method_summary.csv"), row.names = FALSE)
write.csv(aggregated_pips, file.path(out_dir, "combined_aggregated_pips.csv"), row.names = FALSE)
write.csv(agreement_summary, file.path(out_dir, "pip_agreement_with_gagnon.csv"), row.names = FALSE)
write.csv(timing, file.path(out_dir, "paired_timing_with_gagnon.csv"), row.names = FALSE)
writeLines(c(
  paste("baseline:", baseline_dir),
  paste("gagnon:", gagnon_dir),
  "Gagnon comparator: conditional posterior Laplace model probabilities,",
  "Barker locally balanced one-flip weights, posterior-mode/Hessian Gaussian",
  "coefficient proposal, and exact joint MH/RJ correction."
), file.path(out_dir, "manifest.txt"))

print(method_summary)
