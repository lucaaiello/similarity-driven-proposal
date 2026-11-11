#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Rcpp))
args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(name,default) {
  z <- grep(paste0("^--",name,"="),args,value=TRUE)
  if(length(z)) sub(paste0("^--",name,"="),"",z[[1]]) else default
}
budget <- as.numeric(arg_value("budget","2")); chains <- as.integer(arg_value("chains","4"))
out_dir <- arg_value("output",file.path("New_simulations","outputs","online_adaptation_exact"))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
Rcpp::sourceCpp(file.path("New_simulations","online_gaussian_kernels.cpp"),
                cacheDir=file.path("New_simulations",".rcpp_cache"),rebuild=FALSE)

set.seed(20260924L); n <- 100L; p <- 10L; rho <- 0.8
X <- matrix(rnorm(n*p),n,p)
for(j in 2:p) X[,j] <- rho*X[,j-1L]+sqrt(1-rho^2)*X[,j]
X <- unclass(scale(X)); beta <- c(1,-1,rep(0,p-2L)); y <- as.numeric(X%*%beta+rnorm(n))
states <- as.matrix(expand.grid(rep(list(0:1),p)))
logp <- apply(states,1,function(g) evaluate_exact_gaussian_log_target(y,X,as.integer(g)))
prob <- exp(logp-max(logp)); prob <- prob/sum(prob)
exact_pips <- colSums(states*prob)

starts <- list(integer(p),rep(1L,p),as.integer(c(1,0,1,0,1,0,1,0,1,0)))
starts <- rep(starts,length.out=chains)
rows <- list(); estimates <- list(); id <- 0L
for(score in c("f","bic")) for(chain in seq_len(chains)) {
  set.seed(930000L+chain+(score=="bic")*100L)
  fit <- run_optimized_similarity_scaling(
    y,X,score,starts[[chain]],budget,burn=500L,
    lambda_initial=0.7,adapt_lambda=TRUE,balance_directions=TRUE,
    adaptation_metric="diversity",target_effective_neighbors=5,
    adaptation_gain=0.5,adaptation_decay=0.6,adaptation_offset=100,
    uniform_mixture=0.02)
  id <- id+1L
  rows[[id]] <- data.frame(score=score,chain=chain,iterations=fit$iterations,
                           acceptance=fit$postburn_acceptance_rate,
                           final_lambda=fit$final_lambda,
                           mean_effective_neighbors=fit$mean_proposal_effective_neighbors,
                           mean_absolute_pip_error=mean(abs(fit$pips-exact_pips)),
                           max_absolute_pip_error=max(abs(fit$pips-exact_pips)))
  estimates[[id]] <- data.frame(score=score,chain=chain,variable=seq_len(p),
                                 estimate=fit$pips,exact=exact_pips)
}
raw <- do.call(rbind,rows); pip <- do.call(rbind,estimates)
summary <- aggregate(raw[c("acceptance","final_lambda","mean_effective_neighbors",
                            "mean_absolute_pip_error","max_absolute_pip_error")],
                     raw["score"],mean)
write.csv(raw,file.path(out_dir,"chain_summary.csv"),row.names=FALSE)
write.csv(pip,file.path(out_dir,"pip_comparison.csv"),row.names=FALSE)
write.csv(summary,file.path(out_dir,"method_summary.csv"),row.names=FALSE)
saveRDS(list(raw=raw,pip=pip,summary=summary,exact_pips=exact_pips),
        file.path(out_dir,"exact_validation_results.rds"))
print(summary,row.names=FALSE)
