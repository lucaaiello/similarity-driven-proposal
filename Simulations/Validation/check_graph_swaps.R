#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Rcpp))
options(stringsAsFactors=FALSE,warn=1)
args<-commandArgs(trailingOnly=TRUE)
arg_value<-function(name,default){z<-grep(paste0("^--",name,"="),args,value=TRUE);if(length(z))sub(paste0("^--",name,"="),"",z[[1]])else default}
states<-as.integer(arg_value("states","200")); tolerance<-as.numeric(arg_value("tolerance","1e-9"))
out_dir<-arg_value("output",file.path("Simulations/Validation/Results/New_Run/Graph_Swaps"))
source(file.path("Simulations/Validation/gaussian_support.R"));revision_initialize(out_dir,require_audit=FALSE)

Rcpp::sourceCpp(file.path("Simulations/sampler_f_score.cpp"),rebuild=TRUE,verbose=FALSE)
f_swap<-get("audit_f_swap_proposal_cpp",.GlobalEnv)
Rcpp::sourceCpp(file.path("Simulations/sampler_bic_score.cpp"),rebuild=TRUE,verbose=FALSE)
lr_swap<-get("audit_lr_swap_proposal_cpp",.GlobalEnv)

rows<-vector("list",states); set.seed(20269001L)
for(i in seq_len(states)) {
  n<-sample(c(100L,200L),1);p<-sample(c(50L,100L,200L),1)
  X<-scale(matrix(rnorm(n*p),n,p)); y<-as.numeric(.4+X[,1:5]%*%c(.8,-.8,.8,-.8,.8)+rnorm(n))
  gamma<-rbinom(p,1,runif(1,.05,.2)); if(!any(gamma==1))gamma[sample.int(p,1)]<-1
  outgoing<-sample(which(gamma==1),1); candidates<-sample(which(gamma==0),min(8,sum(gamma==0)))
  lambda<-sample(c(.1,.7,1,1.5),1)
  fo<-f_swap(y,cbind(1,X),gamma,outgoing,candidates,lambda,"legacy","informed")
  fn<-f_swap(y,cbind(1,X),gamma,outgoing,candidates,lambda,"rank_one","informed")
  bo<-lr_swap(y,cbind(1,X),gamma,outgoing,candidates,"bic",lambda,"legacy","informed")
  bn<-lr_swap(y,cbind(1,X),gamma,outgoing,candidates,"bic",lambda,"rank_one","informed")
  lo<-lr_swap(y,cbind(1,X),gamma,outgoing,candidates,"lr_power",lambda,"legacy","informed")
  ln<-lr_swap(y,cbind(1,X),gamma,outgoing,candidates,"lr_power",lambda,"rank_one","informed")
  u<-f_swap(y,cbind(1,X),gamma,outgoing,candidates,lambda,"rank_one","uniform")
  rows[[i]]<-data.frame(state=i,n=n,p=p,candidates=length(candidates),lambda=lambda,
    f_max_abs=max(abs(fo-fn)),bic_max_abs=max(abs(bo-bn)),lr_max_abs=max(abs(lo-ln)),
    uniform_max_abs=max(abs(u-rep(1/length(candidates),length(candidates)))))
}
results<-do.call(rbind,rows)
summary<-data.frame(states=nrow(results),tolerance=tolerance,
  f_max_abs=max(results$f_max_abs),bic_max_abs=max(results$bic_max_abs),
  lr_max_abs=max(results$lr_max_abs),uniform_max_abs=max(results$uniform_max_abs))
summary$passed<-max(unlist(summary[c("f_max_abs","bic_max_abs","lr_max_abs","uniform_max_abs")]))<=tolerance
write.csv(results,file.path(out_dir,"statewise_differences.csv"),row.names=FALSE)
write.csv(summary,file.path(out_dir,"summary.csv"),row.names=FALSE)
writeLines(c(paste("states:",states),paste("tolerance:",tolerance),paste("passed:",summary$passed)),file.path(out_dir,"manifest.txt"))
writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"))
print(summary,row.names=FALSE)
if(!summary$passed)stop("Optimized swap probabilities failed the equivalence tolerance")
