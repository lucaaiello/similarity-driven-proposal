source("New_simulations/dm_apt50k_fresh_v1/support.R")
source("New_simulations/dm_apt50k_fresh_v1/runner.R")
apt_environment(); o <- apt_options(c("task", "workers"))
task <- as.integer(apt_get(o, "task", Sys.getenv("SGE_TASK_ID", "0")))
workers <- as.integer(apt_get(o, "workers", Sys.getenv("NSLOTS", "0")))
apt_need(task %in% 1:12 && workers == 10L && workers <= as.integer(Sys.getenv("NSLOTS", "10")),
  "Require task1:12 and10 allocated workers")
z <- long_load(); long_gate()
method <- long_methods[(task - 1L) %% 3L + 1L]; start <- (task - 1L) %/% 3L + 1L
plan <- z$plan; plan$method <- method
plan$seed_base <- plan$seed_base + match(method, long_methods) * 100L
out <- file.path(long_output, method, sprintf("ladder_%02d", start))
long_run(z$inputs, plan, start, out, workers)
cat("FRESH50K TASK COMPLETE", task, method, long_names[start], "\n")
