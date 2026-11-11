# User-run correctness gate. Includes short synthetic smoke chains, NOT real data.
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
source("Simulations/Validation/dm_validation.R")
out <- "Simulations/Validation/Results/New_Run/DM_Coefficients"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(data.frame(passed=FALSE),file.path(out,"summary.csv"),row.names=FALSE)
Rcpp::sourceCpp("Simulations/Dirichlet_Multinomial/dm_sampler.cpp",cacheDir="Simulations/.rcpp_cache")
checks <- list()
check <- function(name, valid, detail="") {
  checks[[length(checks)+1L]] <<- data.frame(check=name,passed=isTRUE(valid),detail=detail)
  if (!isTRUE(valid)) {
    write.csv(do.call(rbind,checks),file.path(out,"checks.csv"),row.names=FALSE)
    stop("Optional-update audit failed: ",name,if (nzchar(detail)) paste0(" (",detail,")") else "")
  }
}
near <- function(x,y,tol=1e-8) max(abs(x-y)) <= tol*(1+max(abs(y)))
set.seed(20261007L)
n <- 12L; P <- 5L; J <- 3L
X <- unclass(scale(matrix(rnorm(n*P),n,P)))
Y <- matrix(sample(0:15,n*J,replace=TRUE),n,J); storage.mode(Y)<-"integer"
alpha <- c(.2,-.3,.4); xb <- matrix(rnorm(n*J,sd=.1),n,J)
for (hyper in list(c(1,9),c(2,3),c(.7,2.5))) for (mask in 0:31) {
  indicator <- as.integer(intToBits(mask))[seq_len(P)]
  expected <- sum(lgamma(indicator+hyper[1])+lgamma(1-indicator+hyper[2]))
  check(paste0("indicator_prior_",paste(hyper,collapse="_"),"_",mask),
    near(audit_dm_indicator_log_prior(indicator,hyper[1],hyper[2]),expected))
}
check("birth_prior_log_odds",near(
  audit_dm_indicator_log_prior(1L)-audit_dm_indicator_log_prior(0L),-log(9)))
check("reproduces_legacy_truncation_difference",near(
  audit_dm_indicator_log_prior(rep(0L,15),legacy_integer_arithmetic=TRUE)-
    audit_dm_indicator_log_prior(rep(0L,15)),-12.027412201221))
check("birth_start_mapping",identical(as.numeric(audit_dm_project_state_start(c(1L,3L),c(.2,-.4),c(1L,2L,3L))),c(.2,0,-.4)))
check("death_start_mapping",near(audit_dm_project_state_start(c(1L,2L,3L),c(.2,.1,-.4),c(1L,3L)),c(.2,-.4)))
check("swap_start_mapping",near(audit_dm_project_state_start(c(1L,3L),c(.2,-.4),c(2L,3L)),c(0,-.4)))
conditional <- function(design,beta,context=xb,budget=3L,warm=TRUE)
  audit_dm_conditional_numerics(design,Y,alpha,context,1L,beta,ncol(design)/n,budget,TRUE,warm)
design <- X[,1:2,drop=FALSE]; b0 <- c(.3,-.1)
first <- conditional(design,b0)
invisible(conditional(design,c(-.4,.5)))
again <- conditional(design,b0)
check("state_start_cache_history_independence",near(first$mode,again$mode))
changed_context <- xb; changed_context[,1] <- rnorm(n)
check("conditional_excludes_fitted_category",near(first$mode,conditional(design,b0,changed_context)$mode))
check("zero_start_ignores_supplied_coefficients",near(conditional(design,b0,warm=FALSE)$mode,
  conditional(design,c(-.4,.5),warm=FALSE)$mode))
# Independent full-kernel likelihood and Gaussian density reference on moderate inputs.
ll <- function(eta) {
  concentration <- exp(eta); sums <- rowSums(concentration)
  sum(lgamma(sums)-lgamma(sums+rowSums(Y)))+sum(lgamma(Y+concentration)-lgamma(concentration))
}
log_normal <- function(value,mean,Q) {
  U <- chol(Q); z <- drop(U %*% (value-mean))
  -.5*(length(value)*log(2*pi)-2*sum(log(diag(U)))+sum(z*z))
}
current <- c(.2,-.1); candidate <- c(.25,-.13)
xb[,1] <- drop(design%*%current)
next_xb <- xb; next_xb[,1] <- drop(design%*%candidate)
terms <- audit_dm_within_model_terms(design,Y,alpha,xb,1L,current,candidate,10,.25)
reverse <- audit_dm_within_model_terms(design,Y,alpha,next_xb,1L,candidate,current,10,.25)
delta <- ll(sweep(next_xb,2,alpha,"+"))-ll(sweep(xb,2,alpha,"+"))
prior <- sum(dnorm(candidate,0,sqrt(10),log=TRUE))-sum(dnorm(current,0,sqrt(10),log=TRUE))
check("within_likelihood_reference",near(terms$likelihood_delta,delta))
check("within_prior_reference",near(terms$prior_delta,prior))
check("within_forward_density_reference",near(terms$log_forward,log_normal(candidate,current,terms$forward_precision)))
check("within_reverse_density_reference",near(terms$log_reverse,log_normal(current,candidate,reverse$forward_precision)))
check("within_ratio_reference",near(terms$ratio,delta+prior+terms$log_reverse-terms$log_forward))
check("within_pair_reversal",near(terms$ratio,-reverse$ratio))
check("within_pair_detailed_balance",near(terms$log_forward+min(0,terms$ratio),
  delta+prior+terms$log_reverse+min(0,reverse$ratio)))
# Fixed birth/death and swap pairs: warm reverse starts from proposed coefficients.
for (scale in c(.02,.25,.5)) {
  fixed <- audit_dm_within_model_terms(design,Y,alpha,xb,1L,current,candidate,10,scale,TRUE)
  back <- audit_dm_within_model_terms(design,Y,alpha,next_xb,1L,candidate,current,10,scale,TRUE)
  Q <- (crossprod(design)/nrow(design)+diag(.1,ncol(design)))/scale^2
  check(paste0("fixed_precision_",scale),near(fixed$forward_precision,Q))
  check(paste0("fixed_density_",scale),near(fixed$log_forward,log_normal(candidate,current,Q)))
  check(paste0("fixed_symmetry_",scale),near(fixed$log_forward,fixed$log_reverse))
  check(paste0("fixed_target_ratio_",scale),near(fixed$ratio,delta+prior))
  check(paste0("fixed_reversal_",scale),near(fixed$ratio,-back$ratio))
  check(paste0("fixed_detailed_balance_",scale),near(fixed$log_forward+min(0,fixed$ratio),
    delta+prior+back$log_forward+min(0,back$ratio)))
  check(paste0("fixed_reference_likelihood_",scale),near(fixed$likelihood_delta,delta))
  check(paste0("fixed_reference_prior_",scale),near(fixed$prior_delta,prior))
}
for (target in list(c(1L,2L,3L),1L,c(2L,3L))) {
  source <- c(1L,3L); saved <- c(.2,-.1); proposed <- rep(.15,length(target))
  forward_start <- audit_dm_project_state_start(source,saved,target)
  reverse_start <- audit_dm_project_state_start(target,proposed,source)
  fm <- conditional(X[,target,drop=FALSE],as.numeric(forward_start))
  rm <- conditional(X[,source,drop=FALSE],as.numeric(reverse_start))
  Qf <- audit_dm_gaussian_numerics(as.numeric(fm$mode),-fm$mode_hessian,proposed)$precision
  Qr <- audit_dm_gaussian_numerics(as.numeric(rm$mode),-rm$mode_hessian,saved)$precision
  check(paste0("warm_forward_density_",paste(target,collapse="_")),near(
    audit_dm_gaussian_numerics(as.numeric(fm$mode),-fm$mode_hessian,proposed)$log_density,
    log_normal(proposed,as.numeric(fm$mode),Qf)))
  check(paste0("warm_reverse_density_",paste(target,collapse="_")),near(
    audit_dm_gaussian_numerics(as.numeric(rm$mode),-rm$mode_hessian,saved)$log_density,
    log_normal(saved,as.numeric(rm$mode),Qr)))
}
neighbors <- lapply(seq_len(P),function(p) setdiff(seq_len(P),p))
initial <- matrix(0L,P,J); initial[1:2,] <- 1L
call_sampler <- function(warm=FALSE,within=FALSE,explicit=TRUE,probability=NULL) {
  settings <- list(X=X,Y=Y,s2=rep(10,J),r2=rep(10,J),a=1,b=9,g_initial=1,c=1,
    runs=60L,burn=30L,thin=1L,adapt_g=TRUE,adaptation_start=5L,window_size=5L,
    adaptation_end=30L,neigh=neighbors,use_move=TRUE,initial_xi=initial,
    coefficient_optimizer_maxit=10L)
  if (explicit) settings <- c(settings,list(coefficient_warm_start=warm,
    within_model_update=within,within_model_attempts=1L,within_model_scale=.25))
  set.seed(20261008L)
  if (!is.null(probability)) settings <- c(settings,list(within_model_fixed_probability=probability,
    within_model_fixed_scale=.25))
  do.call(DM_online_adaptive_VS_RJMCMC,settings)
}
baseline <- call_sampler(explicit=FALSE)
disabled <- call_sampler()
for (field in c("alpha","beta","xi","g_trace","log_posterior_trace"))
  check(paste0("defaults_equal_explicit_off_",field),identical(baseline[[field]],disabled[[field]]))
legacy_within <- call_sampler(TRUE,TRUE)
explicit_zero <- call_sampler(TRUE,TRUE,probability=0)
for (field in c("alpha","beta","xi","g_trace","log_posterior_trace"))
  check(paste0("mixture_zero_equivalence_",field),identical(legacy_within[[field]],explicit_zero[[field]]))
cases <- rbind(expand.grid(warm=c(FALSE,TRUE),within=c(FALSE,TRUE),probability=0),
  data.frame(warm=TRUE,within=TRUE,probability=c(.5,1)))
for (case in seq_len(nrow(cases))) {
  warm <- cases$warm[case]; within <- cases$within[case]; probability <- cases$probability[case]
  fit <- if (!warm&&!within) disabled else call_sampler(warm,within,probability=probability)
  dm_validate_fit(fit)
  check(paste0("mixture_flags_",case),fit$within_model_fixed_probability==probability)
  if (probability==.5) check("both_components_exercised",
    sum(fit$within_fixed_proposed_by_category)>0 && sum(fit$within_curvature_proposed_by_category)>0)
  check(paste0("option_flags_",warm,"_",within),
    identical(fit$coefficient_warm_start,warm)&&identical(fit$within_model_update,within))
  check(paste0("within_attempts_",warm,"_",within),
    if (within) sum(fit$within_model_proposed_by_category)>0 else
      sum(fit$within_model_proposed_by_category)==0)
  for (draw in seq_len(dim(fit$xi)[3])) {
    B <- matrix(0,P,J)
    for (j in seq_len(J)) B[which(fit$xi[,j,draw]==1L),j] <- as.numeric(fit$beta[[draw]][[j]])
    reference <- ll(sweep(X%*%B,2,fit$alpha[draw,],"+"))+
      # Production stores the alpha prior kernel without its fixed normalizing
      # constant. Beta normalizers DO remain because model dimension varies.
      -.5*sum(fit$alpha[draw,]^2/10)+
      sum(dnorm(unlist(fit$beta[[draw]]),0,sqrt(10),log=TRUE))+
      sum(lgamma(fit$xi[,,draw]+1)+lgamma(1-fit$xi[,,draw]+9))
    check(paste0("retained_target_",case,"_",warm,"_",within,"_",draw),near(fit$log_posterior_trace[draw],reference,1e-9),
      sprintf("recorded=%.12g reference=%.12g difference=%.12g",fit$log_posterior_trace[draw],reference,
        fit$log_posterior_trace[draw]-reference))
  }
}
write.csv(do.call(rbind,checks),file.path(out,"checks.csv"),row.names=FALSE)
hashes <- tools::md5sum(dm_source_files)
write.csv(data.frame(file=names(hashes),md5=unname(hashes)),file.path(out,"source_md5.csv"),row.names=FALSE)
write.csv(data.frame(checks=length(checks),passed=TRUE),file.path(out,"summary.csv"),row.names=FALSE)
cat("Optional coefficient-update audit passed:",length(checks),"checks\n")
