# MMV (2013), arXiv1205.1076 equations13-16. Applied here to the LIKELIHOOD
# power path with proper, untempered priors and our trans-dimensional kernels.
# This is their ladder adaptation, not their complete fixed-dimension RW method.
apt_beta <- function(rho) {
  apt_need(all(is.finite(rho)), "Invalid log-gap parameters")
  b <- c(1, exp(-cumsum(exp(rho))))
  apt_need(all(is.finite(b)) && min(b) > 0 && all(diff(b) < 0), "Degenerate temperature ladder")
  b
}
apt_swap_prob <- function(ll, beta) {
  apt_need(length(ll) == length(beta) && all(is.finite(ll)), "Invalid swap likelihoods")
  exp(pmin(0, (head(beta, -1) - tail(beta, -1)) * diff(ll)))
}
apt_rho_update <- function(rho, ll, iteration, target = .234, gain = 1, decay = .6,
                           lower = -12, upper = log(2)) {
  beta <- apt_beta(rho); prob <- apt_swap_prob(ll, beta)
  gamma <- gain * (iteration + 1)^(-decay)
  raw <- rho + gamma * (prob - target)
  list(rho = pmax(lower, pmin(upper, raw)), probabilities = prob,
       projected = raw < lower | raw > upper, gain = gamma)
}
# Preserve the existing category-wise diversity adaptation separately from the
# 23.4% TEMPERATURE-EXCHANGE target. Tuning stays at its temperature slot.
apt_g_update <- function(tuning, effective, plan, informed = TRUE) {
  if (!informed || !plan$adapt_g) return(tuning)
  apt_need(all(is.finite(effective)) && all(effective >= 1 - 1e-10), "Invalid effective-neighbor signal")
  tuning$effective_sum <- tuning$effective_sum + effective
  tuning$effective_count <- tuning$effective_count + 1L
  if (tuning$effective_count == plan$g_window) {
    tuning$g_epoch <- tuning$g_epoch + 1L
    signal <- log(tuning$effective_sum / plan$g_window) - log(plan$g_target)
    eta <- plan$g_scale * max(1, plan$g_offset + tuning$g_epoch - 1)^(-plan$g_decay)
    tuning$g <- exp(pmin(log(10), pmax(log(.05), log(tuning$g) + eta * signal)))
    tuning$effective_sum[] <- 0; tuning$effective_count <- 0L
  }
  tuning
}
apt_flow_update <- function(flow, labels, iteration) {
  L <- length(labels)
  # A completed cold->hot->cold trip for an individual moving physical state.
  for (slot in c(1L, L)) {
    walker <- labels[slot]
    if (slot == 1L && flow$phase[walker] == 0L) {
      flow$phase[walker] <- 1L; flow$departed[walker] <- iteration
    } else if (slot == L && flow$phase[walker] == 1L) flow$phase[walker] <- 2L
    else if (slot == 1L && flow$phase[walker] == 2L) {
      flow$roundtrips[walker] <- flow$roundtrips[walker] + 1L
      flow$durations <- c(flow$durations, iteration - flow$departed[walker])
      flow$phase[walker] <- 1L; flow$departed[walker] <- iteration
    }
  }
  flow
}
