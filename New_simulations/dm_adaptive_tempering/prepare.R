# Metadata and provenance only; no compiler, optimizer, MCMC or network access.
source("New_simulations/dm_adaptive_tempering/support.R")
apt_environment()
x <- readRDS(file.path(apt_root, "inputs/origins.rds"))
original <- "New_simulations/dm_joint_block_pilot"
copied <- c("DM_joint_block.cpp", "joint_kernel.h", "dm_stable_numerics.h")
baseline_matches <- identical(unname(tools::md5sum(file.path(original, copied))),
  unname(tools::md5sum(file.path(apt_root, "baseline", copied))))
input_matches <- identical(unname(tools::md5sum(file.path(original, "inputs/prepared.rds"))),
  unname(tools::md5sum(file.path(apt_root, "inputs/origins.rds"))))
apt_need(baseline_matches && input_matches, "Archived joint pilot/input changed")
plan <- list(schema = "dm_mmv_likelihood_tempering_v1", replicas = 10L, ladders = 4L,
  iterations = 5000L, burn = 2000L, thin = 5L, local_sweeps = 1L,
  target_swap = .234, initial_tmax = 32, gain = 1, decay = .6,
  rho_lower = -12, rho_upper = log(2), adapt_until = 5000L,
  method = "online_similarity_local", joint_scale = 2,
  adapt_g = TRUE, g_window = 25L, g_scale = .5, g_decay = .6, g_offset = 100,
  g_target = min(ncol(x$X), max(5, .025 * ncol(x$X))),
  seed_base = 2026111000L, checkpoint_every = 250L, workers = 4L,
  controls = TRUE, slot_tuning_follows_state = FALSE)
z <- list(plan = plan, source_hashes = apt_hashes(),
  origin_hash = tools::md5sum(file.path(apt_root, "inputs/origins.rds")),
  original_joint_sources = tools::md5sum(file.path(original, copied)))
dest <- file.path(apt_root, "inputs/plan.rds")
if (file.exists(dest)) {
  apt_need(identical(z, readRDS(dest)), "Preserve existing APT plan; sources/settings differ")
} else saveRDS(z, dest)
cat("APT INPUTS READY:10 replicas per ladder;4 independent ladders;4 matched untempered controls.\n")
cat("5000 sweeps/replica, first2000 excluded;MMV temperature adaptation online, target0.234.\n")
cat("Existing diversity-based g adaptation stays separate;original sources/results unchanged.\n")
