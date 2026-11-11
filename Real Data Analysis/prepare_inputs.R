# Pure-R preparation verification only: no fit, sampler or compiler.
source("Real Data Analysis/Sampler/analysis_support.R")
original<-apt_load();supplied<-long_load()
reconstructed<-long_fresh(original$inputs)
stopifnot(identical(reconstructed,supplied$inputs))
cat("PREPARED INPUTS VERIFIED: X/Y, graph and four fresh starts match the supplied inputs.\n")
cat("Original plan and input files were not overwritten.\n")
