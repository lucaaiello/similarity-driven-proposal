#!/usr/bin/env Rscript

# Statewise equality audit for legacy versus rank-one production proposals.
suppressPackageStartupMessages(library(Rcpp))
options(stringsAsFactors=FALSE,warn=1)
args<-commandArgs(trailingOnly=TRUE)
arg_value<-function(name,default){z<-grep(paste0("^--",name,"="),args,value=TRUE);if(length(z))sub(paste0("^--",name,"="),"",z[[1]])else default}
states_per_case<-as.integer(arg_value("states","100"))
tolerance<-as.numeric(arg_value("tolerance","1e-9"))
out_dir<-arg_value("output",file.path("Simulations/Validation/Results/New_Run/Gaussian_Updates"))
source(file.path("Simulations/Validation/gaussian_support.R"));revision_initialize(out_dir,require_audit=FALSE)

compile_exports<-function(path,audit_name){
  Rcpp::sourceCpp(path,rebuild=TRUE,verbose=FALSE)
  get(audit_name,.GlobalEnv)
}
f_audit<-compile_exports(file.path("Simulations/sampler_f_score.cpp"),"audit_f_proposal_cpp")
lr_audit<-compile_exports(file.path("Simulations/sampler_bic_score.cpp"),"audit_lr_proposal_cpp")
Rcpp::sourceCpp(file.path("Simulations/gaussian_kernels.cpp"),rebuild=TRUE,verbose=FALSE)
benchmark_audit<-get("audit_fast_similarity_probabilities",.GlobalEnv)
comparator_audit<-get("audit_fast_comparator_probabilities",.GlobalEnv)

cases<-expand.grid(n=c(50L,100L,200L),p=c(8L,10L,12L),replicate=1:3)
rows<-list();id<-0L
for(case in seq_len(nrow(cases))){
  n<-cases$n[case];p<-cases$p[case];set.seed(20267000L+case)
  rho<-.65;X<-matrix(rnorm(n*p),n,p)
  if(p>1L)for(j in 2:p)X[,j]<-rho*X[,j-1L]+sqrt(1-rho^2)*X[,j]
  X<-unclass(scale(X));y<-as.numeric(.4+X[,unique(round(seq(2,p-1,length.out=3)))]%*%c(1,-1,1)+rnorm(n))
  design<-cbind(1,X)
  for(state in seq_len(states_per_case)){
    gamma<-rbinom(p,1,runif(1,.05,.8))
    exact_lb_difference<-comparator_audit(y,X,gamma,"exact_locally_balanced")$max_absolute_difference
    bic_lb_difference<-comparator_audit(y,X,gamma,"bic_locally_balanced")$max_absolute_difference
    for(balanced in c(FALSE,TRUE))for(lambda in c(.1,.7,1.5)){
      f_old<-f_audit(y,design,gamma,lambda,balanced,if(balanced).02 else 0,"legacy")
      f_new<-f_audit(y,design,gamma,lambda,balanced,if(balanced).02 else 0,"rank_one")
      b_old<-lr_audit(y,design,gamma,"bic",lambda,balanced,if(balanced).02 else 0,"legacy")
      b_new<-lr_audit(y,design,gamma,"bic",lambda,balanced,if(balanced).02 else 0,"rank_one")
      f_benchmark<-benchmark_audit(y,X,gamma,"f",lambda,balanced,if(balanced).02 else 0)$fast
      b_benchmark<-benchmark_audit(y,X,gamma,"bic",lambda,balanced,if(balanced).02 else 0)$fast
      id<-id+1L;rows[[id]]<-data.frame(n=n,p=p,replicate=cases$replicate[case],state=state,
        balanced=balanced,lambda=lambda,f_max_abs=max(abs(f_old-f_new)),
        bic_max_abs=max(abs(b_old-b_new)),
        f_benchmark_max_abs=max(abs(f_new-f_benchmark)),
        bic_benchmark_max_abs=max(abs(b_new-b_benchmark)),
        exact_lb_max_abs=exact_lb_difference,bic_lb_max_abs=bic_lb_difference)
    }
  }
}
results<-do.call(rbind,rows)
write.csv(results,file.path(out_dir,"statewise_differences.csv"),row.names=FALSE)
summary<-data.frame(states=nrow(results),tolerance=tolerance,
  f_max_abs=max(results$f_max_abs),bic_max_abs=max(results$bic_max_abs),
  f_benchmark_max_abs=max(results$f_benchmark_max_abs),
  bic_benchmark_max_abs=max(results$bic_benchmark_max_abs),
  exact_lb_max_abs=max(results$exact_lb_max_abs),bic_lb_max_abs=max(results$bic_lb_max_abs),
  passed=max(results$f_max_abs,results$bic_max_abs,results$f_benchmark_max_abs,
             results$bic_benchmark_max_abs,results$exact_lb_max_abs,
             results$bic_lb_max_abs)<=tolerance)
write.csv(summary,file.path(out_dir,"summary.csv"),row.names=FALSE)
writeLines(c(paste("states_per_case:",states_per_case),paste("cases:",nrow(cases)),
  paste("comparisons:",nrow(results)),paste("tolerance:",tolerance),
  paste("passed:",summary$passed)),file.path(out_dir,"manifest.txt"))
writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"))
print(summary,row.names=FALSE)
if(!summary$passed)stop("Rank-one proposal probabilities failed the equivalence tolerance")
