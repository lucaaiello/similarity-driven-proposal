#!/usr/bin/env Rscript

# Time-budgeted high-dimensional comparison.  This is a scaling and practical
# feasibility experiment; it does not use exact posterior enumeration.
suppressPackageStartupMessages(library(Rcpp))
options(stringsAsFactors=FALSE,warn=1)

args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(name,default) {
  z <- grep(paste0("^--",name,"="),args,value=TRUE)
  if(length(z)) sub(paste0("^--",name,"="),"",z[[1]]) else default
}
mode <- arg_value("mode","quick")
if(!mode %in% c("quick","full")) stop("--mode must be quick or full")
dims <- as.integer(strsplit(arg_value("dimensions",if(mode=="quick") "50,100" else "50,100,250,500,1000"),",",fixed=TRUE)[[1]])
n <- as.integer(arg_value("n","200"))
datasets <- as.integer(arg_value("datasets",if(mode=="quick") "1" else "20"))
chains <- as.integer(arg_value("chains",if(mode=="quick") "2" else "3"))
budget <- as.numeric(arg_value("budget",if(mode=="quick") "0.5" else "5"))
burn <- as.integer(arg_value("burn",if(mode=="quick") "25" else "100"))
max_iterations <- as.integer(arg_value("max-iterations","1000000"))
rho <- as.numeric(arg_value("rho","0.5"))
signals <- as.integer(arg_value("signals","10"))
signal_size <- as.numeric(arg_value("signal-size","0.8"))
out_dir <- arg_value("output","Simulations/Gaussian_Scaling/Results/New_Run")
audit_dir <- arg_value("audit-dir",file.path("Simulations/Validation/Results/New_Run/Gaussian_Updates"))
if(any(!is.finite(dims))||any(dims<signals)||n<10L||datasets<1L||chains<1L||budget<=0||burn<0L)
  stop("Invalid dimensions, n, datasets, chains, budget, burn, or signals")
source(file.path("Simulations/Validation/gaussian_support.R"));revision_initialize(out_dir,require_audit=TRUE)
audit_file<-file.path(audit_dir,"summary.csv")
if(!file.exists(audit_file))stop("Run Simulations/Validation/check_gaussian_updates.R before the scaling experiment")
revision_require_audit(audit_dir)
audit<-read.csv(audit_file)
if(nrow(audit)!=1L||!isTRUE(audit$passed[[1]]))stop("Production/benchmark rank-one equivalence audit did not pass")

Rcpp::sourceCpp(file.path("Simulations/gaussian_kernels.cpp"),
                cacheDir=file.path("Simulations/.rcpp_cache"),rebuild=FALSE)

methods <- data.frame(
  method=c("uniform","f_fixed","bic_fixed","f_online","bic_online",
           "locally_balanced","bic_locally_balanced"),
  engine=c("comparator","production_equivalent_rank_one","production_equivalent_rank_one",
           "production_equivalent_rank_one","production_equivalent_rank_one",
           "comparator","comparator"),
  score=c("uniform","f","bic","f","bic","exact_locally_balanced","bic_locally_balanced"),
  online=c(FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE),stringsAsFactors=FALSE)

ess_basic <- function(x) {
  if(length(x)<20L||!is.finite(var(x))||var(x)==0) return(NA_real_)
  z <- as.numeric(acf(x,lag.max=min(length(x)-1L,1000L),plot=FALSE)$acf)[-1L]
  if(length(z)<2L) return(length(x))
  last <- length(z)-length(z)%%2L
  pairs <- z[seq(1L,last,2L)]+z[seq(2L,last,2L)]
  stop_at <- match(TRUE,pairs<=0,nomatch=length(pairs)+1L)
  kept <- if(stop_at==1L) numeric() else pairs[seq_len(stop_at-1L)]
  min(length(x),length(x)/(1+2*sum(kept)))
}
simulate_data <- function(p,seed) {
  set.seed(seed); X <- matrix(rnorm(n*p),n,p)
  if(p>1L) for(j in 2:p) X[,j] <- rho*X[,j-1L]+sqrt(1-rho^2)*X[,j]
  X <- unclass(scale(X)); active <- unique(as.integer(round(seq(1,p,length.out=signals))))
  beta <- numeric(p); beta[active] <- rep(c(signal_size,-signal_size),length.out=length(active))
  list(X=X,y=as.numeric(0.4+X%*%beta+rnorm(n)),active=active)
}
make_starts <- function(p,active,number,seed) {
  set.seed(seed); z <- replicate(number,integer(p),simplify=FALSE)
  if(number>=2L) z[[2]][active] <- 1L
  if(number>=3L) z[[3]][sample.int(p,length(active))] <- 1L
  if(number>3L) for(i in 4:number) z[[i]][sample.int(p,length(active))] <- 1L
  z
}

rows<-list();pip_rows<-list();id<-0L
for(p in dims) for(dataset in seq_len(datasets)) {
  dat <- simulate_data(p,20264000L+p*100L+dataset)
  initials <- make_starts(p,dat$active,chains,20265000L+p*100L+dataset)
  target_neff <- max(5,0.025*p)
  for(mi in seq_len(nrow(methods))) for(chain in seq_len(chains)) {
    spec <- methods[mi,]; set.seed(20266000L+p*10000L+dataset*100L+mi*10L+chain)
    if(spec$engine=="production_equivalent_rank_one") {
      fit <- run_optimized_similarity_scaling(
        dat$y,dat$X,spec$score,initials[[chain]],budget,burn,max_iterations,
        lambda_initial=0.7,adapt_lambda=spec$online,
        balance_directions=TRUE,adaptation_metric="diversity",
        target_effective_neighbors=target_neff,adaptation_gain=0.5,
        adaptation_decay=0.6,adaptation_offset=100,
        uniform_mixture=0.02)
      acceptance <- fit$postburn_acceptance_rate
      final_lambda <- fit$final_lambda
      mean_neff <- fit$mean_proposal_effective_neighbors
      factorizations <- fit$neighborhood_factorizations
    } else {
      fit <- run_online_gaussian_scaling(
        dat$y,dat$X,spec$score,initials[[chain]],budget,burn,max_iterations,
        lambda=0.7,a=2,d=2,tau=1,a_phi=1,b_phi=3,proposal_engine="rank_one")
      acceptance <- fit$acceptance_rate
      final_lambda <- mean_neff <- NA_real_
      factorizations <- fit$ols_fits
    }
    es <- ess_basic(fit$model_size); et <- ess_basic(fit$log_target)
    top <- order(fit$pips,decreasing=TRUE)[seq_len(min(signals,p))]
    id<-id+1L
    rows[[id]]<-data.frame(
      n=n,dimension=p,dataset=dataset,chain=chain,method=spec$method,
      engine=if(spec$engine=="comparator") "comparator_rank_one" else spec$engine,
      online=spec$online,budget_seconds=budget,total_seconds=fit$total_seconds,
      budget_overrun_seconds=max(0,fit$total_seconds-budget),iterations=fit$iterations,
      retained=fit$retained,iterations_per_second=fit$iterations_per_second,
      accepted_moves_per_second=acceptance*fit$iterations_per_second,
      acceptance_rate=acceptance,ess_model_size=es,ess_log_target=et,
      ess_model_size_per_second=es/fit$total_seconds,
      ess_log_target_per_second=et/fit$total_seconds,
      projected_minutes_per_10000_retained=fit$total_seconds/max(1,fit$retained)*10000/60,
      target_evaluations_per_iteration=fit$target_evaluations/fit$iterations,
      neighborhood_factorizations_per_iteration=factorizations/fit$iterations,
      target_time_share=fit$target_seconds/fit$total_seconds,
      auxiliary_time_share=fit$auxiliary_seconds/fit$total_seconds,
      final_lambda=final_lambda,mean_proposal_effective_neighbors=mean_neff,
      target_effective_neighbors=if(spec$online) target_neff else NA_real_,
      mean_active_pip=mean(fit$pips[dat$active]),
      mean_inactive_pip=mean(fit$pips[-dat$active]),
      top_signal_recall=mean(top%in%dat$active))
    pip_rows[[id]]<-data.frame(n=n,dimension=p,dataset=dataset,chain=chain,
      method=spec$method,variable=seq_len(p),active=seq_len(p)%in%dat$active,pip=fit$pips)
    message(sprintf("p=%d dataset=%d/%d method=%s chain=%d/%d iterations=%d rate=%.1f/s",
                    p,dataset,datasets,spec$method,chain,chains,fit$iterations,fit$iterations_per_second))
  }
}

raw<-do.call(rbind,rows);pips<-do.call(rbind,pip_rows)
numeric_cols<-setdiff(names(raw)[vapply(raw,is.numeric,logical(1))],c("n","dimension","dataset","chain"))
summary_mean<-aggregate(raw[numeric_cols],raw[c("n","dimension","method")],mean,na.rm=TRUE)
summary_se<-aggregate(raw[numeric_cols],raw[c("n","dimension","method")],
  function(x) sd(x,na.rm=TRUE)/sqrt(sum(is.finite(x))))
names(summary_se)[-(1:3)]<-paste0(names(summary_se)[-(1:3)],"_se")
summary<-merge(summary_mean,summary_se,by=c("n","dimension","method"),sort=FALSE)
pip_ranges<-aggregate(pip~n+dimension+dataset+method+variable,pips,function(x) max(x)-min(x))
agreement<-aggregate(pip~n+dimension+method,pip_ranges,mean)
names(agreement)[4]<-"mean_between_chain_pip_range"
summary<-merge(summary,agreement,by=c("n","dimension","method"),sort=FALSE)
exponents<-do.call(rbind,lapply(unique(raw$method),function(m){
  z<-summary[summary$method==m&summary$iterations_per_second>0,]
  data.frame(method=m,empirical_time_exponent=if(nrow(z)>=2L)
    unname(coef(lm(log(1/iterations_per_second)~log(dimension),z))[[2]]) else NA_real_)
}))
write.csv(raw,file.path(out_dir,"chain_metrics.csv"),row.names=FALSE)
write.csv(summary,file.path(out_dir,"method_summary.csv"),row.names=FALSE)
write.csv(pips,file.path(out_dir,"pip_estimates.csv"),row.names=FALSE)
write.csv(exponents,file.path(out_dir,"scaling_exponents.csv"),row.names=FALSE)
writeLines(c(paste("mode:",mode),paste("n:",n),paste("dimensions:",paste(dims,collapse=",")),
  paste("datasets_per_dimension:",datasets),paste("chains:",chains),
  paste("wall_clock_budget_seconds_per_chain_method:",budget),paste("burn_iterations:",burn),
  paste("signals:",signals),paste("rho:",rho),paste("signal_size:",signal_size),
  paste("methods:",paste(methods$method,collapse=",")),
  "target: conjugate Gaussian variable-selection posterior; no exact enumeration",
  "F/BIC engine: rank-one benchmark wrapper, statewise audited against the production rank-one engine",
  "timing: includes target and proposal construction; accepted reverse proposal is cached",
  "stopping: budget checked after each transition and at least burn+50 optimized or burn+10 comparator iterations"),
  file.path(out_dir,"manifest.txt"))
writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"))
saveRDS(list(config=list(mode=mode,n=n,dims=dims,datasets=datasets,chains=chains,budget=budget,burn=burn),
  raw=raw,summary=summary,pips=pips,exponents=exponents),file.path(out_dir,"chain_results.rds"),compress="xz")
print(summary[c("dimension","method","iterations_per_second","acceptance_rate",
  "ess_model_size_per_second","ess_log_target_per_second","projected_minutes_per_10000_retained",
  "top_signal_recall","mean_between_chain_pip_range")],row.names=FALSE)
