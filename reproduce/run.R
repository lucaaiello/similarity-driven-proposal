# Saved-results reproduction only: no MCMC, C++ compilation or cluster access.
args<-commandArgs(TRUE)
getarg<-function(k,default) {
  x<-grep(paste0("^--",k,"="),args,value=TRUE)
  if(length(x))sub(paste0("^--",k,"="),"",x[1])else default
}
unknown<-args[!grepl("^--(output|recompute-parameter-diagnostics)=",args)]
if(length(unknown))stop("Unknown option: ",paste(unknown,collapse=" "))
if(!file.exists("reproduce/run.R"))stop("Run this command from the repository root.")
out<-getarg("output","reproduced")
if(file.exists(out))stop("Choose a NEW output directory; existing results are preserved: ",out)
needed<-c("ggplot2","posterior")
missing<-needed[!vapply(needed,requireNamespace,TRUE,quietly=TRUE)]
if(length(missing))stop("Install required packages first: ",paste(missing,collapse=", "))
summary<-"New_simulations/outputs/dm_apt50k_fresh_slurm_v1/summary/summary.rds"
if(!file.exists(summary)||file.info(summary)$size<1000)
  stop("Real-data saved result missing or an LFS pointer. Run git lfs pull, or extract the complete archive.")
source("reproduce/check_inputs.R")
dir.create(out,recursive=TRUE);out<-normalizePath(out,winslash="/",mustWork=TRUE)
dir.create(file.path(out,"logs"));writeLines("INCOMPLETE",file.path(out,"status.txt"))
exe<-file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe"else"Rscript")
timings<-list()
step<-function(name,script,arguments=character()) {
  cat("Reproducing ",name,"\n",sep="")
  t<-proc.time()[3]
  status<-system2(exe,c(shQuote(script),vapply(arguments,shQuote,"")),
    stdout=file.path(out,"logs",paste0(name,".log")),stderr=file.path(out,"logs",paste0(name,".log")))
  timings[[length(timings)+1L]]<<-data.frame(step=name,seconds=unname(proc.time()[3]-t),exit_code=status)
  write.csv(do.call(rbind,timings),file.path(out,"timings.csv"),row.names=FALSE)
  if(status!=0L)stop("Postprocessing failed: ",name,". Inspect ",file.path(out,"logs",paste0(name,".log")))
}
step("simulation_summaries","New_simulations/make_manuscript_simulation_figures.R",
  c(paste0("--image-dir=",out,"/supporting_figures"),paste0("--summary-dir=",out,"/tables/simulation")))
step("simulation_figures","New_simulations/make_manuscript_simulation_figures_ggplot.R",
  paste0("--image-dir=",out,"/figures"))
step("application","New_simulations/make_manuscript_real_data_figures.R",paste0("--output-root=",out))
step("dm_dataset_diagnostics","New_simulations/summarize_dm_revision_diagnostics.R",
  paste0("--output=",out,"/tables/dm_dataset_diagnostics"))
parameter_input<-"New_simulations/outputs/dm_fitting_depth_20261008_104824/parameter_diagnostics_20261009_221750/parameter_diagnostics.csv"
recompute<-getarg("recompute-parameter-diagnostics","false")
if(!recompute%in%c("true","false"))stop("Use true or false for --recompute-parameter-diagnostics")
if(recompute=="true") {
  step("parameter_diagnostics","New_simulations/summarize_dm_parameter_diagnostics.R",
    paste0("--output=",out,"/tables/recomputed_parameters"))
  parameter_input<-file.path(out,"tables/recomputed_parameters/parameter_diagnostics.csv")
}
step("all_tables","reproduce/tables.R",c(paste0("--output=",out),paste0("--parameters=",parameter_input)))
step("verify_printed_tables","reproduce/check_tables.R",paste0("--output=",out))
file.copy("JCGS_submission/Images/DM_diagram_2.pdf",file.path(out,"figures/DM_diagram_2.pdf"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
write.csv(data.frame(package=needed,version=vapply(needed,function(p)as.character(packageVersion(p)),"")),
  file.path(out,"package_versions.csv"),row.names=FALSE)
writeLines("COMPLETE",file.path(out,"status.txt"))
cat("COMPLETE: ",out,"\nNo sampling, fitting, compilation or cluster jobs were performed.\n",sep="")
