# Run from the repository root. This script DOES execute MCMC when invoked.
args<-commandArgs(TRUE)
opt<-function(k,d) {x<-grep(paste0("^--",k,"="),args,value=TRUE);if(length(x))sub(paste0("^--",k,"="),"",x[1])else d}
mode<-opt("mode","full");stopifnot(mode%in%c("quick","full"))
out<-opt("output","Simulations/Dirichlet_Multinomial/Fitting_Depth/Results/New_Run")
if(file.exists(out))stop("Choose a new output directory; recorded results are preserved.")
dir.create(out,recursive=TRUE)
exe<-file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe"else"Rscript")
invoke<-function(options) {
 status<-system2(exe,c(shQuote("Simulations/Dirichlet_Multinomial/run_dirichlet_multinomial.R"),vapply(options,shQuote,"")))
 if(status!=0L)stop("Simulation stopped; inspect the preceding error.")
}

settings<-if(mode=="full")c("--dataset-ids=1,2,3,4,5","--chains=3","--runs=5000","--burn=2500","--thin=1")else c("--datasets=1","--chains=3","--runs=300","--burn=150","--thin=1")
labels<-c(`1`="One_Step",`5`="Five_Steps",`10`="Ten_Steps")
for(budget in c(1,5,10))invoke(c(paste0("--mode=",mode),"--stage=core",settings,"--methods=optimized_similarity_fixed", "--save-full-fits=true","--cache-score-fits=true",paste0("--score-optimizer-maxit=",budget),paste0("--output=",out,"/",labels[[as.character(budget)]])))
