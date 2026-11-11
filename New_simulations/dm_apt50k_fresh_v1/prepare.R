# File preparation only: no optimizer, compiler or MCMC.
source("New_simulations/dm_apt50k_fresh_v1/support.R")
apt_environment(); old <- apt_load()
dir.create(file.path(long_root, "inputs"), recursive = TRUE, showWarnings = FALSE)
fresh <- long_fresh(old$inputs)
plan <- old$plan
plan$schema <- "fresh_apt_all_edges_50k_burn25k_unthinned_v1"; plan$iterations <- 50000L
plan$burn <- 25000L; plan$thin <- 1L; plan$workers <- 10L; plan$controls <- FALSE
plan$adapt_until <- 50000L; plan$checkpoint_every <- 500L
plan$milestones <- c(12500L, 25000L, 37500L, 50000L)
plan$seed_base <- 2026120000L
plan$exchange_schedule <- "all_adjacent_once_random_permutation_before_local"
plan$initialization <- "fresh_alpha_beta_zero_g_one_covariance_identity"
plan$replica_initialization <- "all_replicas_share_the_ladder_fresh_start"
plan$alpha_covariance_adaptation <- FALSE
plan$methods <- long_methods; plan$start_names <- long_names
plan$max_concurrent_ladders <- 6L
put <- function(z, path) {
  if (file.exists(path)) apt_need(identical(z, readRDS(path)), paste("Preserve existing different input", path))
  else saveRDS(z, path)
}
put(fresh, file.path(long_root, "inputs/fresh.rds"))
put(list(plan = plan, source_hashes = long_hashes(), original_hashes = apt_hashes(),
  input_hash = tools::md5sum(file.path(long_root, "inputs/fresh.rds"))), file.path(long_root, "inputs/plan.rds"))
cat("READY: 3 samplers x4 fresh starts x10 replicas; 50,000 sweeps, burn25,000, no thinning.\n",
  "Nine sequential adjacent exchange attempts per sweep, in random order; 10 one-thread workers per ladder.\n",
  "All tuning reset; no old chain state used. No sampler, optimizer, compiler or submission executed.\n", sep = "")
