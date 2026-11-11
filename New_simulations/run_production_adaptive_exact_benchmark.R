#!/usr/bin/env Rscript

# Replicated exact-posterior audit of the actual production F and BIC samplers,
# comparing matched fixed proposals with online diversity adaptation.
options(stringsAsFactors=FALSE,warn=1)
args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(name,default) {
  z <- grep(paste0("^--",name,"="),args,value=TRUE)
  if(length(z)) sub(paste0("^--",name,"="),"",z[[1]]) else default
}
mode <- arg_value("mode","quick")
if(!mode %in% c("quick","full")) stop("--mode must be quick or full")
design <- arg_value("design","primary")
if(!design %in% c("primary","grid")) stop("--design must be primary or grid")
method_set <- arg_value("method-set","production")
if(!method_set %in% c("production","all")) stop("--method-set must be production or all")
proposal_engine <- arg_value("proposal-engine","rank_one")
if(!proposal_engine %in% c("legacy","rank_one")) stop("--proposal-engine must be legacy or rank_one")
datasets_default <- if(mode=="quick") if(design=="grid") "1" else "2" else if(design=="grid") "20" else "20"
datasets <- as.integer(arg_value("datasets",datasets_default))
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- if(length(script_arg)) sub("^--file=","",script_arg[[1]]) else "run_production_adaptive_exact_benchmark.R"
script_dir <- normalizePath(dirname(script_path),winslash="/",mustWork=TRUE)
repo_dir <- normalizePath(file.path(script_dir,".."),winslash="/",mustWork=TRUE)
output_label <- if(method_set=="all") {
  if(design=="primary") paste0("production_all_methods_exact_",mode) else paste0("production_all_methods_exact_grid_",mode)
} else if(design=="primary") paste0("production_adaptive_exact_",mode) else paste0("production_adaptive_exact_grid_",mode)
if(proposal_engine=="rank_one") output_label <- sub(paste0("_",mode,"$"),paste0("_rank_one_",mode),output_label)
out_dir <- arg_value("output",file.path(script_dir,"outputs",paste0(output_label,"_matched_v2")))
source(file.path("New_simulations","gaussian_revision_support.R"));revision_initialize(out_dir,require_audit=TRUE)

cfg <- if(mode=="quick") list(n=80L,p=8L,active=c(2L,5L,7L),rho=.55,runs=6000L,burn=1500L,chains=2L) else
  list(n=100L,p=10L,active=c(2L,5L,8L),rho=.65,runs=30000L,burn=5000L,chains=4L)
if(design=="grid") {
  cfg$runs <- if(mode=="quick") 4000L else 20000L
  cfg$burn <- if(mode=="quick") 1000L else 5000L
  cfg$chains <- if(mode=="quick") 2L else 3L
}
cfg <- c(cfg,list(intercept=.4,signal=1,sigma=1,a=2,d=2,tau=1,a_phi=1,b_phi=3,lambda=.7,
                  data_seed=20262000L,chain_seed=20263000L))
scenarios <- if(design=="primary") {
  data.frame(n=cfg$n,p=cfg$p)
} else if(mode=="quick") {
  expand.grid(n=c(50L,100L),p=c(6L,8L))
} else {
  expand.grid(n=c(50L,100L,200L),p=c(8L,10L,12L))
}

log_sum_exp <- function(x){m<-max(x);m+log(sum(exp(x-m)))}
all_models <- function(p) t(vapply(0:(2^p-1L),function(s) as.integer(intToBits(s))[seq_len(p)],integer(p)))
generate_data <- function(seed) {
  set.seed(seed); S <- outer(seq_len(cfg$p),seq_len(cfg$p),function(i,j) cfg$rho^abs(i-j))
  x <- scale(matrix(rnorm(cfg$n*cfg$p),cfg$n,cfg$p)%*%chol(S)); b <- numeric(cfg$p)
  b[cfg$active] <- cfg$signal*rep(c(1,-1),length.out=length(cfg$active))
  list(x=unclass(x),y=as.numeric(cfg$intercept+x%*%b+rnorm(cfg$n,sd=cfg$sigma)))
}
log_marginal <- function(bits,design,y) {
  s <- c(1L,which(bits==1L)+1L); G <- crossprod(design[,s,drop=FALSE]); g <- crossprod(design[,s,drop=FALSE],y)
  k <- length(s); P <- diag(1/cfg$tau,k)+G; R <- chol(P)
  m <- backsolve(R,forwardsolve(t(R),g)); astar <- cfg$a+.5*(sum(y^2)-drop(crossprod(m,g)))
  dstar <- cfg$d+length(y)/2
  -.5*length(y)*(log(pi)+log(2))+.5*(-k*log(cfg$tau)-2*sum(log(diag(R))))+
    cfg$d*log(cfg$a)-dstar*log(astar)+lgamma(dstar)-lgamma(cfg$d)
}
enumerate_target <- function(dat) {
  M <- all_models(cfg$p); X <- cbind(1,dat$x); ll <- apply(M,1,log_marginal,design=X,y=dat$y); k <- rowSums(M)
  lp <- lbeta(cfg$a_phi+k,cfg$b_phi+cfg$p-k)-lbeta(cfg$a_phi,cfg$b_phi)
  z <- ll+lp; z <- z-log_sum_exp(z)
  list(models=M,logposterior=z,posterior=exp(z),pips=colSums(M*exp(z)))
}
ess <- function(x) {
  if(length(x)<20L || !is.finite(var(x)) || var(x)==0) return(NA_real_)
  a <- as.numeric(acf(x,plot=FALSE,lag.max=min(length(x)-1L,1000L))$acf)[-1L]
  if(length(a)<2L) return(length(x)); pairs <- a[seq(1,length(a)-1L,2)]+a[seq(2,length(a),2)]
  stop_at <- match(TRUE,pairs<=0,nomatch=length(pairs)+1L); kept <- if(stop_at==1L) numeric() else pairs[seq_len(stop_at-1L)]
  min(length(x),length(x)/(1+2*sum(kept)))
}
compile_sampler <- function(path){Rcpp::sourceCpp(path,rebuild=TRUE,verbose=FALSE);get("var_sel_distance_cpp",.GlobalEnv)}
if(!requireNamespace("Rcpp",quietly=TRUE)||!requireNamespace("RcppArmadillo",quietly=TRUE)) stop("Rcpp packages required")
cat("Compiling production samplers...\n")
f_fun <- compile_sampler(file.path(repo_dir,"Simulations","sampler_F_stat.cpp"))
bic_fun <- compile_sampler(file.path(repo_dir,"Simulations","sampler_LR.cpp"))
methods <- list(
  f_fixed=list(engine="production",fun=f_fun,score=NULL,online=FALSE,seed_offset=100L),
  bic_fixed=list(engine="production",fun=bic_fun,score="bic",online=FALSE,seed_offset=200L),
  f_online=list(engine="production",fun=f_fun,score=NULL,online=TRUE,seed_offset=300L),
  bic_online=list(engine="production",fun=bic_fun,score="bic",online=TRUE,seed_offset=400L))
if(method_set=="all") {
  cat("Compiling exact-reference comparator kernel...\n")
  Rcpp::sourceCpp(file.path(script_dir,"online_gaussian_kernels.cpp"),rebuild=TRUE,verbose=FALSE)
  comparator_fun <- get("run_exact_comparator_kernel",.GlobalEnv)
  methods <- c(methods,list(
    uniform=list(engine="comparator",fun=comparator_fun,score="uniform",online=FALSE,seed_offset=500L),
    locally_balanced=list(engine="comparator",fun=comparator_fun,score="exact_locally_balanced",online=FALSE,seed_offset=600L),
    bic_locally_balanced=list(engine="comparator",fun=comparator_fun,score="bic_locally_balanced",online=FALSE,seed_offset=700L)))
}
chain_rows<-list();pip_rows<-list();trace_rows<-list();exact_rows<-list();id<-0L;tid<-0L
eid<-0L
for(scenario_id in seq_len(nrow(scenarios))) {
 cfg$n <- scenarios$n[scenario_id]; cfg$p <- scenarios$p[scenario_id]
 cfg$active <- unique(as.integer(round(seq(2,cfg$p-1,length.out=3L))))
 initials <- list(rep(0,cfg$p),rep(1,cfg$p),as.numeric(seq_len(cfg$p)%in%cfg$active),as.numeric(seq_len(cfg$p)%%2==0))
 for(dataset in seq_len(datasets)) {
  dat <- generate_data(cfg$data_seed+scenario_id*1000L+dataset); target <- enumerate_target(dat); X <- cbind(1,dat$x)
  eid<-eid+1L; exact_rows[[eid]] <- data.frame(n=cfg$n,p=cfg$p,dataset=dataset,predictor=seq_len(cfg$p),exact_pip=target$pips)
  for(mi in seq_along(methods)) for(chain in seq_len(cfg$chains)) {
    method <- names(methods)[mi]; spec <- methods[[mi]]
    set.seed(cfg$chain_seed+scenario_id*10000L+dataset*1000L+spec$seed_offset+chain)
    if(spec$engine=="production") {
      call <- list(y=dat$y,X=X,neigh=replicate(cfg$p,numeric(),simplify=FALSE),gamma=initials[[chain]],phi=.25,
                   a=cfg$a,d=cfg$d,tau=cfg$tau,a_phi=cfg$a_phi,b_phi=cfg$b_phi,runs=cfg$runs,burn=cfg$burn,
                   thin=1L,start_adapt=100L,end_adapt=1000L,lambda_bd=cfg$lambda,lambda_move=1,
                   adapt=spec$online,move_step=FALSE,proposal_engine=proposal_engine)
      if(!is.null(spec$score)) call$score_type <- spec$score
      call <- c(call,list(adaptation_mode="diversity",balanced_directions=TRUE,
                                          uniform_mixture=.02,target_effective_fraction=.025,
                                          adaptation_gain=.5,adaptation_decay=.6,adaptation_offset=100))
    } else {
      call <- list(y=dat$y,predictors=dat$x,method=spec$score,gamma=as.integer(initials[[chain]]),
                   runs=cfg$runs,burn=cfg$burn,a=cfg$a,d=cfg$d,tau=cfg$tau,
                   a_phi=cfg$a_phi,b_phi=cfg$b_phi,proposal_engine=proposal_engine)
    }
    cat(sprintf("n=%d p=%d dataset=%d/%d method=%s chain=%d/%d\n",cfg$n,cfg$p,dataset,datasets,method,chain,cfg$chains))
    started <- proc.time()[["elapsed"]]; fit <- do.call(spec$fun,call); elapsed <- proc.time()[["elapsed"]]-started
    G <- fit$gamma; codes <- rowSums(sweep(G,2,2^(0:(cfg$p-1L)),`*`)); counts <- tabulate(codes+1L,nbins=nrow(target$models))
    est <- counts/length(codes); pip <- colMeans(G); size <- rowSums(G); logtarget <- target$logposterior[codes+1L]
    if(spec$engine=="production") {
      acceptance_rate <- fit$acc_rate_bd
      lambda_tail <- tail(fit$lambda_chain,max(1L,floor(length(fit$lambda_chain)/4L)))
      neff_tail <- tail(fit$proposal_effective_neighbors,max(1L,floor(length(fit$proposal_effective_neighbors)/4L)))
      final_lambda <- tail(fit$lambda_chain,1)
      lambda_range <- max(lambda_tail)-min(lambda_tail)
      mean_neff <- mean(neff_tail)
      target_neff <- fit$target_effective_neighbors
    } else {
      acceptance_rate <- fit$acceptance_rate
      final_lambda <- lambda_range <- mean_neff <- target_neff <- NA_real_
    }
    id<-id+1L; chain_rows[[id]]<-data.frame(n=cfg$n,p=cfg$p,dataset=dataset,method=method,chain=chain,online=spec$online,
      engine=if(spec$engine=="production") paste0("production_",proposal_engine) else paste0("comparator_",proposal_engine),
      elapsed_seconds=elapsed,acceptance_rate=acceptance_rate,unique_models=sum(counts>0),
      ess_model_size=ess(size),ess_log_target=ess(logtarget),ess_model_size_per_second=ess(size)/elapsed,
      ess_log_target_per_second=ess(logtarget)/elapsed,mean_abs_pip_error=mean(abs(pip-target$pips)),
      max_abs_pip_error=max(abs(pip-target$pips)),total_variation_distance=.5*sum(abs(est-target$posterior)),
      final_lambda=final_lambda,lambda_last_quarter_range=lambda_range,
      mean_last_quarter_effective_neighbors=mean_neff,target_effective_neighbors=target_neff)
    pip_rows[[id]]<-data.frame(n=cfg$n,p=cfg$p,dataset=dataset,method=method,chain=chain,predictor=seq_len(cfg$p),
                               exact_pip=target$pips,estimated_pip=pip,absolute_error=abs(pip-target$pips))
    if(spec$online){keep<-unique(c(seq(1,cfg$runs,by=100L),cfg$runs));tid<-tid+1L
      trace_rows[[tid]]<-data.frame(n=cfg$n,p=cfg$p,dataset=dataset,method=method,chain=chain,iteration=keep,
                                    lambda=fit$lambda_chain[keep],effective_neighbors=fit$proposal_effective_neighbors[keep])}
  }
 }
}
chains_df<-do.call(rbind,chain_rows);pips_df<-do.call(rbind,pip_rows);exact_df<-do.call(rbind,exact_rows)
traces_df<-if(length(trace_rows))do.call(rbind,trace_rows) else data.frame()
num<-setdiff(names(chains_df)[vapply(chains_df,is.numeric,logical(1))],c("n","p","dataset","chain"))
means<-aggregate(chains_df[num],chains_df[c("n","p","method")],mean,na.rm=TRUE)
ses<-aggregate(chains_df[num],chains_df[c("n","p","method")],function(x)sd(x,na.rm=TRUE)/sqrt(sum(is.finite(x))))
names(ses)[-(1:3)]<-paste0(names(ses)[-(1:3)],"_se")
summary_df<-merge(means,ses,by=c("n","p","method"),sort=FALSE)
ranges<-aggregate(estimated_pip~n+p+dataset+method+predictor,pips_df,function(x)max(x)-min(x))
mr<-aggregate(estimated_pip~n+p+method,ranges,max);names(mr)[4]<-"max_between_chain_pip_range"
summary_df<-merge(summary_df,mr,by=c("n","p","method"),sort=FALSE)
overall_means<-aggregate(chains_df[num],chains_df["method"],mean,na.rm=TRUE)
timing_baseline<-chains_df[chains_df$method=="locally_balanced",
  c("n","p","dataset","chain","elapsed_seconds","ess_model_size_per_second","ess_log_target_per_second")]
names(timing_baseline)[5:7]<-c("locally_balanced_seconds","locally_balanced_ess_model_size_per_second",
  "locally_balanced_ess_log_target_per_second")
timing_comparisons<-merge(chains_df,timing_baseline,by=c("n","p","dataset","chain"),sort=FALSE)
timing_comparisons$runtime_ratio_to_locally_balanced<-timing_comparisons$elapsed_seconds/timing_comparisons$locally_balanced_seconds
timing_comparisons$ess_model_size_per_second_ratio_to_locally_balanced<-
  timing_comparisons$ess_model_size_per_second/timing_comparisons$locally_balanced_ess_model_size_per_second
timing_comparisons$ess_log_target_per_second_ratio_to_locally_balanced<-
  timing_comparisons$ess_log_target_per_second/timing_comparisons$locally_balanced_ess_log_target_per_second
write.csv(chains_df,file.path(out_dir,"chain_summary.csv"),row.names=FALSE)
write.csv(summary_df,file.path(out_dir,"method_summary.csv"),row.names=FALSE)
write.csv(overall_means,file.path(out_dir,"overall_method_summary.csv"),row.names=FALSE)
write.csv(pips_df,file.path(out_dir,"pip_estimates.csv"),row.names=FALSE)
write.csv(exact_df,file.path(out_dir,"exact_pips.csv"),row.names=FALSE)
write.csv(traces_df,file.path(out_dir,"adaptive_trajectories_thinned.csv"),row.names=FALSE)
write.csv(timing_comparisons,file.path(out_dir,"paired_timing_comparisons.csv"),row.names=FALSE)
saveRDS(list(config=cfg,design=design,scenarios=scenarios,datasets=datasets,chains=chains_df,summary=summary_df,
             overall=overall_means,pips=pips_df,traces=traces_df),
        file.path(out_dir,"production_adaptive_exact_results.rds"),compress="xz")
writeLines(c(paste("mode:",mode),paste("design:",design),paste("method_set:",method_set),
             paste("production_proposal_engine:",proposal_engine),paste("datasets_per_scenario:",datasets),
             paste("scenarios:",paste0("n",scenarios$n,"_p",scenarios$p,collapse=",")),paste("runs:",cfg$runs),
             paste("burn:",cfg$burn),paste("chains:",cfg$chains),
             "target: exact production conjugate posterior with phi integrated to Beta-binomial configuration prior",
             paste("methods:",paste(names(methods),collapse=",")),
             "timing: all proposal-score construction is included; accepted reverse proposals are cached as the next forward proposal",
             "paired_timing_baseline: locally_balanced on the same dataset and chain start"),file.path(out_dir,"manifest.txt"))
writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"))
print(summary_df[c("n","p","method","acceptance_rate","ess_model_size","ess_log_target","mean_abs_pip_error",
                   "max_abs_pip_error","total_variation_distance","final_lambda",
                   "mean_last_quarter_effective_neighbors","max_between_chain_pip_range")],row.names=FALSE)
