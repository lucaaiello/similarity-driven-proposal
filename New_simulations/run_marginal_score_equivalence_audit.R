#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Rcpp))
options(stringsAsFactors=FALSE,warn=1)
args<-commandArgs(trailingOnly=TRUE)
arg_value<-function(name,default){z<-grep(paste0("^--",name,"="),args,value=TRUE);if(length(z))sub(paste0("^--",name,"="),"",z[[1]])else default}
states<-as.integer(arg_value("states","200"));tolerance<-as.numeric(arg_value("tolerance","1e-10"))
out_dir<-arg_value("output",file.path("New_simulations","outputs","marginal_score_equivalence"))
source(file.path("New_simulations","gaussian_revision_support.R"));revision_initialize(out_dir,require_audit=FALSE)
Rcpp::sourceCpp(file.path("New_simulations","online_gaussian_kernels.cpp"),cacheDir=file.path("New_simulations",".rcpp_cache"),rebuild=FALSE)

reference_probabilities<-function(y,X,gamma,lambda,balanced,mix){
  n<-length(y);p<-ncol(X);yc<-y-mean(y);rss0<-sum(yc^2)
  evidence<-vapply(seq_len(p),function(j){x<-X[,j]-mean(X[,j]);explained<-sum(x*yc)^2/sum(x^2);residual<-max(rss0-explained,1e-14);f<-max(explained/(residual/max(n-2,1)),0);max(-pf(f,1,max(n-2,1),lower.tail=FALSE,log.p=TRUE)/log(10),0)},numeric(1))
  active_sum<-sum(evidence[gamma==1]);candidate_k<-sum(gamma)+ifelse(gamma==1,-1,1)
  score<-ifelse(candidate_k>0,(active_sum+ifelse(gamma==1,-evidence,evidence))/sqrt(pmax(candidate_k,1)),0)
  logw<-score^lambda;prob<-numeric(p)
  if(!balanced){z<-exp(logw-max(logw));prob<-z/sum(z);if(mix>0)prob<-(1-mix)*prob+mix/p
  }else{birth<-which(gamma==0);death<-which(gamma==1);bm<-if(length(birth)&&length(death)).5 else as.numeric(length(birth)>0);dm<-if(length(birth)&&length(death)).5 else as.numeric(length(death)>0);if(length(birth)){z<-exp(logw[birth]-max(logw[birth]));prob[birth]<-bm*z/sum(z)};if(length(death)){z<-exp(logw[death]-max(logw[death]));prob[death]<-dm*z/sum(z)};prob<-(1-mix)*prob;if(length(birth))prob[birth]<-prob[birth]+bm*mix/length(birth);if(length(death))prob[death]<-prob[death]+dm*mix/length(death)}
  prob
}

set.seed(20272001L);rows<-vector("list",states)
for(i in seq_len(states)){n<-sample(c(50L,100L,200L),1);p<-sample(c(8L,10L,50L,200L),1);X<-scale(matrix(rnorm(n*p),n,p));y<-as.numeric(.4+X[,seq_len(min(3,p))]%*%rep(c(1,-1,1),length.out=min(3,p))+rnorm(n));gamma<-integer(p);model_size<-sample.int(min(p,n-3L,30L)+1L,1L)-1L;if(model_size>0L)gamma[sample.int(p,model_size)]<-1L;lambda<-sample(c(.1,.5,.7,1,1.5),1);balanced<-sample(c(FALSE,TRUE),1);mix<-if(balanced).02 else 0;reference<-reference_probabilities(y,X,gamma,lambda,balanced,mix);optimized<-evaluate_fast_similarity_probabilities(y,X,gamma,"marginal",lambda,balanced,mix);rows[[i]]<-data.frame(state=i,n=n,p=p,model_size=model_size,lambda=lambda,balanced=balanced,max_abs=max(abs(reference-optimized)),sum_error=abs(sum(optimized)-1))}
raw<-do.call(rbind,rows);finite<-is.finite(raw$max_abs)&is.finite(raw$sum_error);summary<-data.frame(states=nrow(raw),valid_finite_states=sum(finite),tolerance=tolerance,max_abs=if(any(finite))max(raw$max_abs[finite])else Inf,max_sum_error=if(any(finite))max(raw$sum_error[finite])else Inf);summary$passed<-all(finite)&&max(summary$max_abs,summary$max_sum_error)<=tolerance
write.csv(raw,file.path(out_dir,"statewise_differences.csv"),row.names=FALSE);write.csv(summary,file.path(out_dir,"summary.csv"),row.names=FALSE);writeLines(c(paste("states:",states),paste("valid_finite_states:",sum(finite)),paste("tolerance:",tolerance),paste("passed:",summary$passed)),file.path(out_dir,"manifest.txt"));writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"));print(summary,row.names=FALSE);if(!isTRUE(summary$passed))stop("Marginal-score optimized probabilities failed the independent reference audit")
