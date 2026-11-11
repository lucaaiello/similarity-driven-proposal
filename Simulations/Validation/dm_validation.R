# Shared hard gate before checkpointing or reporting real-data diagnostics.
dm_numerical_engine <- "stable_rising_factorial_v1"
dm_fit_engine <- "deterministic_zero_start_exact_context_cache_v1"
dm_source_files <- c(
  "Simulations/Dirichlet_Multinomial/dm_sampler.cpp",
  "Simulations/Dirichlet_Multinomial/dm_numerics.h",
  "Simulations/Validation/dm_validation.R"
)
dm_validate_fit <- function(fit, label = "chain") {
  if (!identical(fit$indicator_prior_engine,"floating_log_gamma_v1"))
    stop(label, ": obsolete indicator-prior arithmetic; do not reuse this checkpoint for the corrected target")
  required <- c("alpha", "beta", "xi", "log_posterior_trace", "g_trace")
  if (!all(required %in% names(fit))) stop(label, ": missing retained output")
  if (!length(fit$log_posterior_trace) ||
      any(!is.finite(fit$log_posterior_trace))) {
    stop(label, ": nonfinite/missing posterior trace; no convergence or PIP summary is valid")
  }
  if (any(!is.finite(fit$alpha)) || any(!is.finite(unlist(fit$beta))) ||
      any(!is.finite(fit$xi)) || any(!fit$xi %in% c(0L, 1L)) ||
      any(!is.finite(fit$g_trace))) stop(label, ": invalid retained state")
  dims <- dim(fit$xi)
  if (length(dims) != 3L || dims[3] != length(fit$log_posterior_trace) ||
      length(fit$beta) != dims[3] || !identical(dim(fit$alpha), dims[2:3][2:1])) {
    stop(label, ": inconsistent retained dimensions")
  }
  for (draw in seq_len(dims[3])) {
    if (length(fit$beta[[draw]]) != dims[2]) stop(label, ": beta category mismatch")
    sizes <- vapply(fit$beta[[draw]], length, integer(1))
    if (any(sizes != colSums(fit$xi[, , draw, drop = FALSE], dims = 1L))) {
      stop(label, ": coefficient/indicator size mismatch")
    }
  }
  if (!identical(fit$numerical_engine, dm_numerical_engine) ||
      !fit$fit_engine %in% c(dm_fit_engine, "deterministic_state_start_exact_context_cache_v1")) {
    stop(label, ": obsolete sampler engine; use a fresh output directory and rerun")
  }
  if (identical(fit$fit_engine, "deterministic_state_start_exact_context_cache_v1") &&
      !isTRUE(fit$coefficient_warm_start)) stop(label, ": state-start engine missing warm-start flag")
  if (isTRUE(fit$coefficient_warm_start) && identical(fit$fit_engine, dm_fit_engine))
    stop(label, ": warm-start flag inconsistent with engine")
  if (isTRUE(fit$within_model_update)) {
    proposed <- fit$within_model_proposed_by_category
    accepted <- fit$within_model_accepted_by_category
    if (length(proposed) != dims[2] || length(accepted) != dims[2] ||
        any(!is.finite(c(proposed,accepted))) || any(accepted < 0 | accepted > proposed))
      stop(label, ": invalid within-model counters")
    if (!is.null(fit$within_model_fixed_probability)) {
      probability <- fit$within_model_fixed_probability
      if (length(probability)!=1L || !is.finite(probability) || probability<0 || probability>1 ||
          length(fit$within_model_fixed_scale)!=1L || !is.finite(fit$within_model_fixed_scale) ||
          fit$within_model_fixed_scale<=0) stop(label,": invalid mixture settings")
      fp <- fit$within_fixed_proposed_by_category; fa <- fit$within_fixed_accepted_by_category
      cp <- fit$within_curvature_proposed_by_category; ca <- fit$within_curvature_accepted_by_category
      if (any(lengths(list(fp,fa,cp,ca))!=dims[2]) || any(!is.finite(c(fp,fa,cp,ca))) ||
          any(c(fp,fa,cp,ca)<0) || any(fa>fp | ca>cp) ||
          any(fp+cp!=proposed | fa+ca!=accepted) ||
          (probability==0 && any(fp!=0)) || (probability==1 && any(cp!=0)))
        stop(label,": invalid mixture counters")
    }
  }
  invisible(TRUE)
}
