# Run from the repository root. This script DOES execute MCMC when invoked.
args<-commandArgs(TRUE)
opt<-function(k,d) {x<-grep(paste0("^--",k,"="),args,value=TRUE);if(length(x))sub(paste0("^--",k,"="),"",x[1])else d}
mode<-opt("mode","full");stopifnot(mode%in%c("quick","full"))
out<-opt("output","Simulations/Dirichlet_Multinomial/Comparison/Results/New_Run")
if(file.exists(out))stop("Choose a new output directory; recorded results are preserved.")
dir.create(out,recursive=TRUE)
exe<-file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe"else"Rscript")
invoke<-function(options) {
 status<-system2(exe,c(shQuote("Simulations/Dirichlet_Multinomial/run_dirichlet_multinomial.R"),vapply(options,shQuote,"")))
 if(status!=0L)stop("Simulation stopped; inspect the preceding error.")
}

settings<-if(mode=="full")c("--datasets=20","--chains=3","--runs=5000","--burn=2500","--thin=5")else c("--datasets=1","--chains=2","--runs=300","--burn=150","--thin=1")
for(methods in c("uniform_flip,optimized_similarity_fixed,online_similarity_flip,online_similarity_local","gagnon_laplace_rj")) {
 label<-if(methods=="gagnon_laplace_rj")"Gagnon_Laplace"else"Similarity_Methods"
 invoke(c(paste0("--mode=",mode),"--stage=core",settings,paste0("--methods=",methods),paste0("--output=",out,"/",label)))
}
