#!/usr/bin/env Rscript
# Exhaustive frozen-kernel audit + checks of actual production sampler traces.
options(stringsAsFactors=FALSE,warn=1)
args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(name,default) {
  z <- grep(paste0("^--",name,"="),args,value=TRUE)
  if(length(z)) sub(paste0("^--",name,"="),"",z[[1]]) else default
}
out_dir <- arg_value("output","New_simulations/outputs/matched_gaussian_revision_audit")
source("New_simulations/gaussian_revision_support.R")
revision_initialize(out_dir,FALSE)
suppressPackageStartupMessages(library(Rcpp))
exports <- list()
for(score in c("f","bic")) {
  file <- if(score=="f") "Simulations/sampler_F_stat.cpp" else "Simulations/sampler_LR.cpp"
  sourceCpp(file,rebuild=TRUE)
  exports[[score]] <- list(sampler=get("var_sel_distance_cpp"),
    proposal=get(if(score=="f")"audit_f_proposal_cpp" else "audit_lr_proposal_cpp"),
    target=get(if(score=="f")"audit_f_log_target_cpp" else "audit_lr_log_target_cpp"))
}
sourceCpp("New_simulations/online_gaussian_kernels.cpp",rebuild=TRUE)
reference_target <- function(g,X,y,a=2,d=2,tau=1,ap=1,bp=3) {
  Z <- X[,c(1,which(g==1)+1),drop=FALSE]
  P <- crossprod(Z)+diag(1/tau,ncol(Z)); C <- chol(P); rhs <- crossprod(Z,y)
  mean <- backsolve(C,forwardsolve(t(C),rhs)); astar <- a+.5*(sum(y^2)-sum(mean*rhs))
  ds <- d+length(y)/2; k <- sum(g)
  -length(y)/2*log(2*pi)+.5*(-ncol(Z)*log(tau)-2*sum(log(diag(C))))+
    d*log(a)-ds*log(astar)+lgamma(ds)-lgamma(d)+lbeta(ap+k,bp+length(g)-k)-lbeta(ap,bp)
}
rows <- list()
check <- function(name,error,tolerance=1e-8) {
  passed <- is.finite(error) && error<=tolerance
  rows[[length(rows)+1L]] <<- data.frame(check=name,error=error,tolerance=tolerance,passed=passed)
  message(name,": ",format(error),if(passed)" PASS" else " FAIL")
}
set.seed(20261007)
p <- 5L
states <- t(vapply(0:(2^p-1L),function(i)as.integer(intToBits(i))[1:p],integer(p)))
for(design in c("regular","dense","duplicate")) {
  n <- if(design=="dense")4L else 20L
  x <- matrix(rnorm(n*p),n,p)
  if(design=="duplicate")x[,2] <- x[,1]
  X <- cbind(1,x); y <- .4+x[,1]+rnorm(n)
  target <- apply(states,1,function(g)reference_target(g,X,y))
  posterior <- exp(target-max(target));posterior <- posterior/sum(posterior)
  proposals <- list()
  for(score in c("f","bic","marginal","uniform","exact_locally_balanced","bic_locally_balanced")) {
    if(score%in%c("f","bic")) {
      fun <- exports[[score]]$proposal
      Q <- t(apply(states,1,function(g) {
        call <- list(y=y,X=X,gamma=g,lambda=.7,balanced_directions=TRUE,
          uniform_mixture=.02,proposal_engine="rank_one")
        if(score=="bic")call$score_type <- "bic"
        do.call(fun,call)
      }))
      cpp_target <- apply(states,1,function(g)exports[[score]]$target(y,X,g))
      check(paste(design,score,"production target"),max(abs(target-cpp_target)))
      Qfast <- t(apply(states,1,function(g)evaluate_fast_similarity_probabilities(
        y,x,as.integer(g),score,.7,TRUE,.02)))
      check(paste(design,score,"production/benchmark"),max(abs(Q-Qfast)))
    } else if(score=="marginal") {
      Q <- t(apply(states,1,function(g)evaluate_fast_similarity_probabilities(
        y,x,as.integer(g),score,.7,TRUE,.02)))
      cpp_target <- apply(states,1,function(g)evaluate_exact_gaussian_log_target(y,x,as.integer(g)))
    } else {
      Q <- t(apply(states,1,function(g)audit_fast_comparator_probabilities(
        y,x,as.integer(g),score)$fast))
      cpp_target <- apply(states,1,function(g)evaluate_exact_gaussian_log_target(y,x,as.integer(g)))
    }
    check(paste(design,score,"positive full support"),as.numeric(any(!is.finite(Q)|Q<=0)),0)
    check(paste(design,score,"normalization"),max(abs(rowSums(Q)-1)))
    K <- matrix(0,nrow(states),nrow(states))
    for(i in seq_len(nrow(states)))for(j in seq_len(p)) {
      next_i <- bitwXor(i-1L,bitwShiftL(1L,j-1L))+1L
      accept <- exp(min(0,cpp_target[next_i]-cpp_target[i]+log(Q[next_i,j])-log(Q[i,j])))
      K[i,next_i] <- Q[i,j]*accept
    }
    diag(K) <- 1-rowSums(K)
    flux <- K*posterior
    check(paste(design,score,"detailed balance"),max(abs(flux-t(flux))))
    check(paste(design,score,"stationarity"),max(abs(drop(posterior%*%K)-posterior)))
    if(score%in%c("f","bic","marginal")) {
      for(i in seq_len(nrow(states))) {
        g <- states[i,]
        if(sum(g)>0 && sum(g)<p)check(paste(design,score,i,"birth mass"),abs(sum(Q[i,g==0])-.5))
      }
    }
    saveRDS(list(states=states,target=target,posterior=posterior,Q=Q,K=K),
      file.path(out_dir,paste0(design,"_",score,"_transition.rds")))
  }
}
# Actual sampler: no update before start_adapt => exact fixed/adaptive equality.
set.seed(123);x <- matrix(rnorm(40*p),40,p);X <- cbind(1,x);y <- .4+x[,1]+rnorm(40)
for(score in c("f","bic")) {
  call <- list(y=y,X=X,neigh=replicate(p,numeric(),simplify=FALSE),gamma=numeric(p),
    phi=.25,a=2,d=2,tau=1,a_phi=1,b_phi=3,runs=600L,burn=100L,thin=1L,
    start_adapt=10000L,lambda_bd=.7,adapt=FALSE,move_step=FALSE,
    adaptation_mode="diversity",balanced_directions=TRUE,uniform_mixture=.02,
    proposal_engine="rank_one")
  if(score=="bic")call$score_type <- "bic"
  set.seed(321);fixed <- do.call(exports[[score]]$sampler,call)
  call$adapt <- TRUE
  set.seed(321);frozen <- do.call(exports[[score]]$sampler,call)
  check(paste(score,"fixed/frozen adaptive path"),max(abs(fixed$gamma-frozen$gamma)),0)
  call$start_adapt <- 100L
  set.seed(321);online <- do.call(exports[[score]]$sampler,call)
  check(paste(score,"adaptation actually updates"),as.numeric(all(online$lambda_chain==.7)),0)
  for(regime in c("fixed","online")) {
    fit <- get(regime)
    expected <- apply(fit$gamma,1,function(g)reference_target(g,X,y))
    check(paste(score,regime,"recorded collapsed log target"),max(abs(fit$loglik-expected)))
    qcall <- list(y=y,X=X,gamma=as.numeric(fit$final_gamma),
      lambda=tail(fit$lambda_chain,1),balanced_directions=TRUE,
      uniform_mixture=.02,proposal_engine="rank_one")
    if(score=="bic")qcall$score_type <- "bic"
    fresh <- do.call(exports[[score]]$proposal,qcall)
    check(paste(score,regime,"final proposal cache"),max(abs(fresh-fit$final_proposal)))
  }
  # Auxiliary phi is not part of the collapsed model acceptance ratio.
  call$adapt <- FALSE;call$phi <- .9
  set.seed(321);other_phi <- do.call(exports[[score]]$sampler,call)
  check(paste(score,"initial phi independence"),max(abs(fixed$gamma-other_phi$gamma)),0)
}
summary <- do.call(rbind,rows)
write.csv(summary,file.path(out_dir,"summary.csv"),row.names=FALSE)
if(!all(summary$passed))stop("Matched Gaussian revision audit failed; do not run benchmarks.")
message("All matched Gaussian revision checks passed.")
