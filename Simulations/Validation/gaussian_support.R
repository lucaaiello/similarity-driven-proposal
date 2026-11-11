# Run provenance and fail-fast checks shared by corrected Gaussian experiments.
revision_sources <- function() {
  paths <- c("Simulations/sampler_f_score.cpp", "Simulations/sampler_bic_score.cpp",
    "Simulations/gaussian_kernels.cpp",
    "Simulations/gaussian_model.h",
    "Simulations/Validation/gaussian_support.R",
    "Simulations/Validation/check_gaussian_kernels.R",
    "Simulations/Validation/check_gaussian_updates.R",
    "Simulations/Validation/check_graph_swaps.R",
    "Simulations/Validation/check_marginal_scores.R")
  data.frame(path=paths, md5=unname(tools::md5sum(paths)))
}
revision_record <- function(out_dir) {
  write.csv(revision_sources(), file.path(out_dir,"source_hashes.csv"),row.names=FALSE)
  script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  saveRDS(list(time=Sys.time(), command=commandArgs(FALSE),
    script_hash=tools::md5sum(script), sources=revision_sources(),
    session=sessionInfo(), threads=Sys.getenv(c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"))),
    file.path(out_dir,"provenance.rds"))
}
revision_require_audit <- function(directory) {
  summary_path <- file.path(directory,"summary.csv")
  hash_path <- file.path(directory,"source_hashes.csv")
  if (!file.exists(summary_path) || !file.exists(hash_path))
    stop("Missing current revision audit: ",directory)
  result <- read.csv(summary_path)
  if (!nrow(result) || !all(result$passed %in% TRUE))
    stop("Revision audit has not passed: ",directory)
  previous <- read.csv(hash_path)
  current <- revision_sources()
  if (!identical(previous$path,current$path) || !identical(previous$md5,current$md5))
    stop("Sources changed since audit; rerun audits: ",directory)
}
revision_initialize <- function(out_dir, require_audit=TRUE) {
  if (dir.exists(out_dir) && length(list.files(out_dir,all.files=TRUE,no..=TRUE)))
    stop("Refusing to overwrite existing results: ",out_dir)
  if (require_audit) revision_require_audit(arg_value("revision-audit",
    file.path("Simulations/Validation/Results/New_Run/Gaussian_Kernels")))
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  revision_record(out_dir)
}
