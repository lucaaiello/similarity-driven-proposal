# Post-process saved full fits only: this script never launches the sampler.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = "") {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
root <- "Simulations/Dirichlet_Multinomial/Fitting_Depth/Results"
inputs <- strsplit(option("inputs", paste(file.path(root, c("One_Step","Five_Steps","Ten_Steps")), collapse = ";")), ";", fixed = TRUE)[[1]]
output <- option("output")
inspect <- "--inspect" %in% args
self_test <- "--self-test" %in% args
max_groups <- as.integer(option("max-groups", "0"))
if (is.na(max_groups) || max_groups < 0L) stop("Invalid --max-groups")
if (nzchar(output) && dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) stop("Use a fresh output directory")

reconstruct_beta <- function(fit) {
  dims <- dim(fit$xi)
  result <- array(0, dims)
  if (length(fit$beta) != dims[3]) stop("Coefficient/indicator draw count mismatch")
  for (draw in seq_len(dims[3])) {
    if (length(fit$beta[[draw]]) != dims[2]) stop("Coefficient category count mismatch")
    for (j in seq_len(dims[2])) {
      active <- which(fit$xi[, j, draw] == 1L)
      b <- as.numeric(fit$beta[[draw]][[j]])
      if (length(b) != length(active)) stop("Active coefficient/indicator mismatch")
      result[active, j, draw] <- b
    }
  }
  result
}

if (!requireNamespace("posterior", quietly = TRUE) && !inspect) stop("Install the R package 'posterior' before running diagnostics")
diagnose <- function(x) {
  stopifnot(is.matrix(x), ncol(x) >= 2L, nrow(x) >= 20L, all(is.finite(x)))
  constant <- apply(x, 2, function(z) all(z == z[1]))
  switches <- apply(x, 2, function(z) sum(diff(z) != 0))
  # A constant chain supplies no evidence of agreement, even if all chains match.
  informative <- !any(constant)
  result <- c(rhat = NA_real_, bulk_ess = NA_real_, tail_ess = NA_real_,
              mean_ess = NA_real_, mcse_mean = NA_real_)
  if (informative) {
    result <- c(rhat = posterior::rhat(x), bulk_ess = posterior::ess_bulk(x),
                tail_ess = posterior::ess_tail(x), mean_ess = posterior::ess_mean(x),
                mcse_mean = posterior::mcse_mean(x))
  }
  c(result, informative = informative, constant_chains = sum(constant),
    min_switches = min(switches), mean = mean(x), chain_mean_range = diff(range(colMeans(x))))
}
if (self_test) {
  set.seed(123)
  good <- diagnose(matrix(rnorm(3000), 1000, 3))
  flat <- diagnose(matrix(0, 100, 3))
  partial <- matrix(rnorm(300), 100, 3); partial[, 1] <- 0
  stopifnot(is.finite(good["rhat"]), is.na(flat["rhat"]),
            flat["constant_chains"] == 3, is.na(diagnose(partial)["rhat"]))
  fake <- list(xi = array(c(1, 0, 0, 1), c(2, 1, 2)), beta = list(list(2), list(3)))
  stopifnot(identical(as.numeric(reconstruct_beta(fake)), c(2, 0, 0, 3)))
  cat("Self-test passed; no saved chains were changed.\n")
  quit(status = 0)
}

rows <- list(); provenance <- list(); group_count <- 0L
for (input in inputs) {
  full_dir <- file.path(input, "full_fits")
  files <- sort(list.files(full_dir, pattern = "\\.rds$", full.names = TRUE))
  if (!length(files)) stop("No full fits in ", full_dir, "; aggregate CSVs cannot reconstruct parameter ESS or Rhat")
  if (inspect) {
    z <- readRDS(files[1])
    cat("\nInput:", input, "\nFull fits:", length(files), "\n")
    print(z$settings); print(dim(z$fit$xi)); print(dim(z$fit$alpha))
    print(z$fit[c("numerical_engine", "fit_engine", "indicator_prior_engine")])
    next
  }
  saved <- lapply(files, readRDS)
  keys <- vapply(saved, function(z) paste(z$dataset, z$method, z$quality, sep = "/"), character(1))
  groups <- split(seq_along(saved), keys)
  for (indices in groups) {
    if (max_groups && group_count >= max_groups) break
    records <- saved[indices]
    records <- records[order(vapply(records, function(z) z$chain, numeric(1)))]
    if (length(records) < 3L || anyDuplicated(vapply(records, function(z) z$chain, numeric(1)))) stop("Require at least three distinct chains per dataset")
    dims <- dim(records[[1]]$fit$xi)
    if (length(dims) != 3L || dims[3] < 20L) stop("Invalid retained indicator dimensions")
    for (z in records) {
      if (!identical(dim(z$fit$xi), dims) ||
          !identical(dim(z$fit$alpha), c(dims[3], dims[2])) ||
          any(!is.finite(z$fit$xi)) || any(!z$fit$xi %in% c(0, 1)) ||
          any(!is.finite(z$fit$alpha)) ||
          !identical(z$data$xi, records[[1]]$data$xi) ||
          !identical(z$data$beta, records[[1]]$data$beta)) stop("Invalid or unmatched full fits")
      provenance[[length(provenance) + 1L]] <- data.frame(input = input, dataset = z$dataset,
        method = z$method, chain = z$chain, score_steps = z$settings$score_optimizer_maxit,
        retained_draws = dims[3], burn = z$settings$burn, thin = z$settings$thin,
        numerical_engine = if (is.null(z$fit$numerical_engine)) "unrecorded" else z$fit$numerical_engine,
        fit_engine = if (is.null(z$fit$fit_engine)) "unrecorded" else z$fit$fit_engine,
        indicator_prior_engine = if (is.null(z$fit$indicator_prior_engine)) "unrecorded" else z$fit$indicator_prior_engine)
    }
    beta <- lapply(records, function(z) reconstruct_beta(z$fit))
    if (any(!is.finite(unlist(beta)))) stop("Nonfinite saved coefficients")
    add <- function(x, quantity, predictor = NA_integer_, category = NA_integer_, active = NA) {
      rows[[length(rows) + 1L]] <<- cbind(data.frame(input = input,
        dataset = records[[1]]$dataset, method = records[[1]]$method,
        score_steps = records[[1]]$settings$score_optimizer_maxit,
        quantity = quantity, predictor = predictor, category = category,
        generating_active = active, chains = length(records), draws_per_chain = dims[3]),
        as.data.frame(as.list(diagnose(x))))
    }
    for (j in seq_len(dims[2])) {
      add(vapply(records, function(z) z$fit$alpha[, j], numeric(dims[3])), "intercept", category = j)
      for (p in seq_len(dims[1])) {
        active <- as.logical(records[[1]]$data$xi[p, j])
        add(vapply(records, function(z) as.numeric(z$fit$xi[p, j, ]), numeric(dims[3])), "indicator", p, j, active)
        # Inactive coefficients are zero: do not condition away excluded draws.
        add(vapply(beta, function(b) as.numeric(b[p, j, ]), numeric(dims[3])), "coefficient", p, j, active)
      }
    }
    group_count <- group_count + 1L
    cat("Completed", basename(input), "dataset", records[[1]]$dataset, "\n")
  }
  if (max_groups && group_count >= max_groups) break
}
if (inspect) {
  cat("\nposterior package available:", requireNamespace("posterior", quietly = TRUE), "\n")
  quit(status = 0)
}
parameters <- do.call(rbind, rows)
parameters$parameter_class <- ifelse(parameters$quantity == "intercept", "intercept",
  paste(parameters$quantity, ifelse(parameters$generating_active, "active", "inactive"), sep = "_"))
finite_summary <- function(x, fun) if (any(is.finite(x))) fun(x[is.finite(x)]) else NA_real_
keys <- interaction(parameters$input, parameters$dataset, parameters$parameter_class, drop = TRUE)
dataset_summary <- do.call(rbind, lapply(split(parameters, keys), function(z) data.frame(
  input = z$input[1], dataset = z$dataset[1], score_steps = z$score_steps[1], parameter_class = z$parameter_class[1],
  parameters = nrow(z), informative = sum(z$informative == 1),
  any_constant_chain = sum(z$constant_chains > 0),
  rhat_above_1_01 = sum(z$rhat > 1.01, na.rm = TRUE),
  rhat_above_1_05 = sum(z$rhat > 1.05, na.rm = TRUE),
  max_rhat = finite_summary(z$rhat, max),
  min_bulk_ess = finite_summary(z$bulk_ess, min), min_tail_ess = finite_summary(z$tail_ess, min),
  max_mcse_mean = finite_summary(z$mcse_mean, max), max_chain_mean_range = max(z$chain_mean_range))))
flags <- parameters[parameters$constant_chains > 0 | !is.finite(parameters$rhat) |
                     parameters$rhat > 1.01, ]
print(dataset_summary, row.names = FALSE)
cat("\nConstant chains are reported as uninformative, not as converged.\n")
cat("MCSE is for the posterior mean; indicator MCSE is therefore PIP MCSE.\n")
cat("These diagnostics do not validate the target or establish convergence.\n")
if (nzchar(output)) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  write.csv(parameters, file.path(output, "parameter_diagnostics.csv"), row.names = FALSE)
  write.csv(dataset_summary, file.path(output, "dataset_parameter_summary.csv"), row.names = FALSE)
  write.csv(flags, file.path(output, "flagged_parameters.csv"), row.names = FALSE)
  write.csv(do.call(rbind, provenance), file.path(output, "provenance.csv"), row.names = FALSE)
  selected <- parameters[parameters$quantity != "coefficient" | parameters$generating_active, ]
  table_keys <- interaction(selected$score_steps, selected$quantity, drop = TRUE)
  table_summary <- do.call(rbind, lapply(split(selected, table_keys), function(z) data.frame(
    score_steps = z$score_steps[1], quantity = z$quantity[1],
    informative = sum(is.finite(z$rhat)), total = nrow(z),
    rhat_above_1_01 = sum(z$rhat > 1.01, na.rm = TRUE),
    max_rhat = finite_summary(z$rhat, max),
    min_bulk_ess = finite_summary(z$bulk_ess, min),
    max_mcse_mean = finite_summary(z$mcse_mean, max))))
  table_summary <- table_summary[order(table_summary$score_steps, table_summary$quantity), ]
  fmt <- function(x, digits) if (is.finite(x)) formatC(x, format = "f", digits = digits) else "--"
  labels <- c(coefficient = "Active coefficients", indicator = "Inclusion indicators", intercept = "Intercepts")
  latex_rows <- vapply(seq_len(nrow(table_summary)), function(i) {
    z <- table_summary[i, ]
    paste0(paste(z$score_steps, labels[z$quantity], paste0(z$informative, "/", z$total),
      fmt(z$max_rhat, 3), fmt(z$min_bulk_ess, 0), fmt(z$max_mcse_mean, 4), sep = " & "), " \\\\")
  }, character(1))
  write.csv(table_summary, file.path(output, "manuscript_parameter_summary.csv"), row.names = FALSE)
  writeLines(latex_rows, file.path(output, "table_rows.tex"))
  used <- unlist(lapply(inputs, function(d) list.files(file.path(d, "full_fits"), pattern = "\\.rds$", full.names = TRUE)))
  hashes <- tools::md5sum(used)
  write.csv(data.frame(file = names(hashes), md5 = unname(hashes)), file.path(output, "input_md5.csv"), row.names = FALSE)
  cat("\nWritten to", output, "\n")
}
