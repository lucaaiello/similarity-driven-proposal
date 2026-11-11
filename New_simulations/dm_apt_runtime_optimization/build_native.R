# Deterministic source construction only. No compiler, optimizer or MCMC.
# Prefix numerical helpers and the statistical sweep are copied verbatim from
# the immutable reference; only initialization/output scaffolding is replaced.
root <- "New_simulations/dm_apt_runtime_optimization"
reference <- "New_simulations/dm_adaptive_tempering/DM_tempered.cpp"
dir.create(file.path(root, "baseline"), recursive = TRUE, showWarnings = FALSE)
for (name in c("DM_tempered.cpp", "joint_kernel.h", "dm_stable_numerics.h")) {
  from <- file.path("New_simulations/dm_adaptive_tempering", name)
  to <- file.path(root, "baseline", name)
  if (file.exists(to)) stopifnot(identical(unname(tools::md5sum(from)), unname(tools::md5sum(to))))
  else stopifnot(file.copy(from, to))
}
for (name in c("joint_kernel.h", "dm_stable_numerics.h")) {
  from <- file.path("New_simulations/dm_adaptive_tempering", name); to <- file.path(root, name)
  if (file.exists(to)) stopifnot(identical(unname(tools::md5sum(from)), unname(tools::md5sum(to))))
  else stopifnot(file.copy(from, to))
}
lines <- readLines(reference, warn = FALSE, encoding = "UTF-8")
loc <- function(pattern) { z <- grep(pattern, lines, fixed = TRUE); stopifnot(length(z) == 1L); z }
entry <- loc("Rcpp::List DM_tempered_VS_RJMCMC(")
stopifnot(lines[entry - 1L] == "// [[Rcpp::export]]")
prefix <- lines[seq_len(entry - 2L)]
begin <- loc('  if (!Sigma_alpha.is_finite()) Rcpp::stop("Nonfinite alpha covariance at iteration %d",iter_id);')
end <- loc("  // After the loop: increment counter if any beta changed")
body <- lines[seq.int(begin, end - 1L)]
# This compact API reproduces the reference worker regime: burn=0/adapt_g=false
# inside the native sweep, with diversity and temperatures updated by R.
a <- grep("  // Freeze the intercept proposal as well as g for retained sampling.", body, fixed = TRUE)
b <- grep("  // (beta, xi) update //", body, fixed = TRUE)
stopifnot(length(a) == 1L, length(b) == 1L, b > a)
# Retain the section divider following the dead burn-in branch.
body <- body[-seq.int(a, b - 2L)]
tail_start <- grep("    bool in_adaptation_window = adapt_g &&", body, fixed = TRUE)
stopifnot(length(tail_start) == 1L)
body <- c(body[seq_len(tail_start - 1L)], "  }")

prologue <- r"---(

// Immutable data/graph plus reusable physical-state views. A context contains
// no RNG and no private adaptive tuning: those remain in the R slot state.
// Reconstruct after serialization from the same immutable inputs.
class APTPersistentContext {
 public:
  arma::mat X;
  arma::imat Y;
  int n,J,P;
  arma::vec s2,r2;
  double a=1.0,b=9.0,c=1.0;
  bool informed=true,use_move=true,balance_directions=true,cache_score_fits=true;
  bool coefficient_warm_start=true,within_model_update=true,joint_update=true;
  int n_move_attempts=1,within_model_attempts=1,joint_attempts=1,laplace_optimizer_maxit=25;
  double lambda_move=1.0,uniform_mixture=.02;
  double within_model_scale=.25,within_model_fixed_probability=.5,within_model_fixed_scale=.05;
  std::string discrete_proposal="similarity";
  std::vector<std::vector<int>> neigh_index;
  arma::ivec y_sum_i,PP;
  arma::vec alpha,log_prior_beta,log_prior_xi,lambda,gamma_sum_i;
  arma::imat xi;
  std::vector<std::string> config;
  std::vector<arma::uvec> selected;
  std::vector<arma::mat> X_list;
  std::vector<arma::vec> beta;
  std::vector<std::unordered_map<std::string,arma::vec>> beta_hat_memory;
  arma::mat Xbeta,alpha_mat,eta,gamma,prop_unif;
  double log_lik=0,log_prior_alpha=0;
  std::uint64_t steps=0;

  APTPersistentContext(const arma::mat& data_X,const arma::imat& data_Y,
      const Rcpp::List& neighbors,const std::string& method):X(data_X),Y(data_Y),
      n(Y.n_rows),J(Y.n_cols),P(X.n_cols),s2(J,arma::fill::value(10.0)),
      r2(J,arma::fill::value(10.0)),PP(J),alpha(J),log_prior_beta(J),log_prior_xi(J),
      lambda(J),gamma_sum_i(n),xi(P,J),config(J),selected(J),X_list(J),beta(J),
      beta_hat_memory(J),Xbeta(n,J),alpha_mat(n,J),eta(n,J),gamma(n,J),
      prop_unif(P,J,arma::fill::ones) {
    if(method!="online_similarity_local"&&method!="online_similarity_flip"&&method!="uniform_flip")
      Rcpp::stop("Persistent context supports the three current similarity/uniform methods only");
    informed=method!="uniform_flip";use_move=method=="online_similarity_local";
    balance_directions=informed;uniform_mixture=informed ? .02 : 0.0;
    if(n<1||J<1||P<1||X.n_rows!=Y.n_rows||!X.is_finite()||Y.min()<0)
      Rcpp::stop("Invalid persistent DM data");
    for(int i=0;i<n;++i) if(arma::accu(arma::conv_to<arma::rowvec>::from(Y.row(i)))>
        std::numeric_limits<int>::max()) Rcpp::stop("DM row count exceeds integer range");
    if(use_move) {
      if(neighbors.size()!=static_cast<R_xlen_t>(P))Rcpp::stop("Invalid persistent neighborhood size");
      neigh_index.resize(P);
      for(int p=0;p<P;++p) {
        arma::ivec neigh_p=Rcpp::as<arma::ivec>(neighbors[p]);
        for(arma::uword idx=0;idx<neigh_p.n_elem;++idx) {
          int q=neigh_p(idx)-1;
          if(q>=0&&q<P&&q!=p)neigh_index[p].push_back(q);
        }
      }
      for(int p=0;p<P;++p) {
        std::vector<int> copy=neigh_index[p];
        for(int q:copy)if(q>=0&&q<P&&q!=p)neigh_index[q].push_back(p);
      }
      for(int p=0;p<P;++p) {
        std::vector<int>& v=neigh_index[p];std::sort(v.begin(),v.end());
        v.erase(std::unique(v.begin(),v.end()),v.end());
      }
    }
    y_sum_i=arma::sum(Y,1);
    for(int j=0;j<J;++j)prop_unif.col(j)/=P;
  }

  Rcpp::List step(const Rcpp::List& physical,const arma::vec& g,
      const arma::mat& proposal_covariance,double inverse_temperature,
      int coefficient_optimizer_maxit,int score_optimizer_maxit,double joint_scale) {
    auto begin_step=std::chrono::steady_clock::now();
    if(!std::isfinite(inverse_temperature)||inverse_temperature<=0||inverse_temperature>1||
       coefficient_optimizer_maxit<1||score_optimizer_maxit<1||!std::isfinite(joint_scale)||joint_scale<=0)
      Rcpp::stop("Invalid persistent step settings");
    dm_numeric_diagnostics::rejected_ratios=0;
    dm_numeric_diagnostics::precision_projections=0;
    dm_numeric_diagnostics::optimizer_fallbacks=0;
    alpha=Rcpp::as<arma::vec>(physical["alpha"]);
    xi=Rcpp::as<arma::imat>(physical["xi"]);
    beta=Rcpp::as<std::vector<arma::vec>>(physical["beta"]);
    arma::mat Sigma_alpha=proposal_covariance;
    if(alpha.n_elem!=static_cast<arma::uword>(J)||xi.n_rows!=static_cast<arma::uword>(P)||
       xi.n_cols!=static_cast<arma::uword>(J)||beta.size()!=static_cast<size_t>(J)||
       g.n_elem!=static_cast<arma::uword>(J)||!alpha.is_finite()||!g.is_finite()||g.min()<=0||
       Sigma_alpha.n_rows!=static_cast<arma::uword>(J)||Sigma_alpha.n_cols!=static_cast<arma::uword>(J)||
       !Sigma_alpha.is_finite()||arma::any(arma::vectorise(xi)<0)||arma::any(arma::vectorise(xi)>1))
      Rcpp::stop("Invalid persistent physical/tuning state");
    arma::mat check_chol;
    if(!arma::chol(check_chol,arma::symmatu(Sigma_alpha)))Rcpp::stop("Persistent alpha covariance is not SPD");
    PP=arma::sum(xi,0).t();
    // Rebase exactly as the old pilot_start path does on EVERY invocation.
    // This deliberately preserves its floating-point refresh semantics.
    for(int j=0;j<J;++j) {
      selected[j]=arma::find(xi.col(j)>0);X_list[j]=X.cols(selected[j]);
      config[j]=vectorToKey(arma::conv_to<arma::ivec>::from(xi.col(j)));
      if(beta[j].n_elem!=selected[j].n_elem||!beta[j].is_finite())Rcpp::stop("Persistent beta/indicator mismatch");
      Xbeta.col(j)=X_list[j]*beta[j];
      log_prior_beta(j)=-(PP(j)/2.0)*(std::log(2.0*M_PI)+std::log(r2(j)))-.5*arma::dot(beta[j],beta[j])/r2(j);
      log_prior_xi(j)=dm_indicator_log_prior(xi.col(j),a,b);
      lambda(j)=c*static_cast<double>(PP(j))/static_cast<double>(n);
      beta_hat_memory[j].clear();if(PP(j)>0)beta_hat_memory[j][config[j]]=beta[j];
    }
    alpha_mat.each_row()=alpha.t();eta=alpha_mat+Xbeta;gamma=arma::trunc_exp(eta);
    gamma_sum_i=arma::sum(gamma,1);log_lik=dm_stable_log_kernel(Y,eta);
    log_prior_alpha=-.5*arma::accu(arma::square(alpha)/s2);
    if(!std::isfinite(log_lik)||!std::isfinite(log_prior_alpha)||
       !log_prior_beta.is_finite()||!log_prior_xi.is_finite()||!Xbeta.is_finite())
      Rcpp::stop("Nonfinite initial persistent DM state");
    Rcpp::checkUserInterrupt();
    const int iter_id=1;
    int accept_alpha=0,accept_betaxi_global=0,accept_move_global=0,prop_move_global=0;
    arma::ivec joint_proposed(J,arma::fill::zeros),joint_accepted(J,arma::fill::zeros);
    arma::ivec within_proposed(J,arma::fill::zeros),within_accepted(J,arma::fill::zeros);
    arma::ivec within_fixed_proposed(J,arma::fill::zeros),within_fixed_accepted(J,arma::fill::zeros);
    arma::ivec within_curvature_proposed(J,arma::fill::zeros),within_curvature_accepted(J,arma::fill::zeros);
    arma::ivec accept_betaxi(J,arma::fill::zeros),accept_move(J,arma::fill::zeros),prop_move(J,arma::fill::zeros);
    arma::vec proposal_effective_this_iter(J,arma::fill::zeros);
    arma::vec proposal_entropy_this_iter(J,arma::fill::zeros),proposal_max_this_iter(J,arma::fill::zeros);
)---"

epilogue <- r"---(
    // Same full end-of-sweep target refresh as the old runs=1/thin=1 worker.
    log_lik=dm_stable_log_kernel(Y,alpha_mat+Xbeta);
    if(!std::isfinite(log_lik))Rcpp::stop("Nonfinite retained persistent DM likelihood");
    const double recorded_target=log_lik+log_prior_alpha+arma::accu(log_prior_beta)+arma::accu(log_prior_xi);
    if(!std::isfinite(recorded_target))Rcpp::stop("Nonfinite retained persistent DM posterior");
    ++steps;
    long double sb=0,sx=0;
    for(int j=0;j<J;++j){sb+=log_prior_beta(j);sx+=log_prior_xi(j);}
    const double lp=log_prior_alpha+static_cast<double>(sb)+static_cast<double>(sx);
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin_step).count();
    Rcpp::NumericVector counts=Rcpp::NumericVector::create(
      Rcpp::Named("alpha")=static_cast<double>(accept_alpha),
      Rcpp::Named("flip")=static_cast<double>(arma::accu(accept_betaxi)),
      Rcpp::Named("local")=static_cast<double>(arma::accu(accept_move)),
      Rcpp::Named("joint")=static_cast<double>(arma::accu(joint_accepted)),
      Rcpp::Named("projections")=static_cast<double>(dm_numeric_diagnostics::precision_projections.load()));
    return Rcpp::List::create(
      Rcpp::Named("physical")=Rcpp::List::create(Rcpp::Named("alpha")=alpha,
        Rcpp::Named("xi")=xi,Rcpp::Named("beta")=Rcpp::wrap(beta)),
      Rcpp::Named("ll")=log_lik,Rcpp::Named("lp")=lp,
      Rcpp::Named("effective")=proposal_effective_this_iter,
      Rcpp::Named("counts")=counts,Rcpp::Named("local_seconds")=seconds,
      Rcpp::Named("nonfinite_mh_ratios_rejected")=dm_numeric_diagnostics::rejected_ratios.load(),
      Rcpp::Named("optimizer_zero_fallbacks")=dm_numeric_diagnostics::optimizer_fallbacks.load(),
      Rcpp::Named("prior_alpha")=log_prior_alpha,Rcpp::Named("prior_beta")=log_prior_beta,
      Rcpp::Named("prior_xi")=log_prior_xi);
  }
};

// [[Rcpp::export]]
SEXP apt_context_create(const arma::mat& X,const arma::imat& Y,
    const Rcpp::List& neighbors,std::string method="online_similarity_local") {
  Rcpp::XPtr<APTPersistentContext> ptr(new APTPersistentContext(X,Y,neighbors,method),true);
  R_SetExternalPtrTag(ptr,Rf_install("APTPersistentContext_v1"));
  return ptr;
}

// Sigma_alpha is the ALREADY temperature-scaled proposal covariance supplied by
// the old apt_args/pilot_start wrapper; this API does not rescale it a second time.
// [[Rcpp::export]]
Rcpp::List apt_context_step(SEXP context,const Rcpp::List& physical,
    const arma::vec& g,const arma::mat& Sigma_alpha,double inverse_temperature,
    int coefficient_optimizer_maxit=25,int score_optimizer_maxit=1,double joint_scale=2.0) {
  if(TYPEOF(context)!=EXTPTRSXP||R_ExternalPtrAddr(context)==nullptr||
     R_ExternalPtrTag(context)!=Rf_install("APTPersistentContext_v1"))
    Rcpp::stop("Invalid or serialized native context; recreate it from immutable inputs");
  Rcpp::XPtr<APTPersistentContext> ptr(context);
  return ptr->step(physical,g,Sigma_alpha,inverse_temperature,
    coefficient_optimizer_maxit,score_optimizer_maxit,joint_scale);
}

// No RNG, optimizer or sampler; supports ownership/restart regression checks.
// [[Rcpp::export]]
Rcpp::List apt_context_info(SEXP context) {
  if(TYPEOF(context)!=EXTPTRSXP||R_ExternalPtrAddr(context)==nullptr||
     R_ExternalPtrTag(context)!=Rf_install("APTPersistentContext_v1"))
    Rcpp::stop("Invalid or serialized native context; recreate it from immutable inputs");
  Rcpp::XPtr<APTPersistentContext> ptr(context);
  return Rcpp::List::create(Rcpp::Named("n")=ptr->n,Rcpp::Named("P")=ptr->P,
    Rcpp::Named("J")=ptr->J,Rcpp::Named("steps")=static_cast<double>(ptr->steps),
    Rcpp::Named("informed")=ptr->informed,Rcpp::Named("use_move")=ptr->use_move,
    Rcpp::Named("engine")="persistent_immutable_context_compact_step_v1");
}
)---"

destination <- file.path(root, "DM_persistent.cpp")
if (file.exists(destination)) stop("Preserve existing generated DM_persistent.cpp before rebuilding")
writeLines(c("// GENERATED computational copy; original sources/results remain immutable.",
  paste0("// Reference MD5: ", unname(tools::md5sum(reference))), prefix,
  prologue, body, epilogue), destination, useBytes = TRUE)
saveRDS(list(reference = reference, reference_md5 = tools::md5sum(reference),
  native_md5 = tools::md5sum(destination), copied_sweep_lines = length(body),
  removed_branches = c("dead native alpha warmup", "dead native g/post-run diagnostic accumulation")),
  file.path(root, "native_generation.rds"))
cat("NATIVE SOURCE GENERATED; no compiler, optimizer or sampler executed.\n")
