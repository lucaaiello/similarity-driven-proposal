#ifndef DM_OPTIONAL_JOINT_KERNEL_H
#define DM_OPTIONAL_JOINT_KERNEL_H
// Included after the existing conditional-Hessian and Gaussian helpers.
// Fixed-model, full-rank MH: theta=(alpha_j,beta_j), no optimizer or Jacobian.
namespace dm_joint {
struct Terms {
  double likelihood_delta, alpha_prior_delta, beta_prior_delta;
  double log_forward, log_reverse, ratio;
  arma::mat raw_forward, raw_reverse, forward_precision, reverse_precision;
};
inline arma::vec coordinates(double alpha,const arma::vec& beta) {
  arma::vec result(beta.n_elem+1);result(0)=alpha;
  if(beta.n_elem)result.tail(beta.n_elem)=beta;
  return result;
}
inline arma::mat raw_precision(const arma::mat& design,const arma::imat& Y,
    const arma::vec& alpha,const arma::ivec& pp,const arma::mat& xb,int j,
    const arma::vec& theta,double s2,double r2,double inverse_temperature=1.0) {
  arma::mat Z=arma::join_rows(arma::ones(design.n_rows,1),design);
  arma::vec zero_alpha=alpha;zero_alpha(j)=0;
  arma::ivec augmented_pp=pp;augmented_pp(j)=Z.n_cols;
  arma::mat raw=-dm_reg_hessian(Z,Y,zero_alpha.t(),augmented_pp,xb,j,theta,1.0/r2);
  // The augmented helper uses 1/r2 for every coordinate; correct alpha's
  // actual prior precision independently. This is NOT the k/n PMLE penalty.
  raw(0,0)+=1.0/s2-1.0/r2;
  if(inverse_temperature!=1.0) {
    arma::mat prior=arma::eye(raw.n_rows,raw.n_cols)/r2;prior(0,0)=1.0/s2;
    raw=inverse_temperature*(raw-prior)+prior;
  }
  return arma::symmatu(raw);
}
inline arma::mat precision(const arma::mat& raw,double scale) {
  return dm_gaussian_precision(raw)*(static_cast<double>(raw.n_rows)/(scale*scale));
}
inline Terms terms(const arma::mat& design,const arma::imat& Y,
    const arma::vec& alpha,const arma::ivec& pp,const arma::mat& xb,int j,
    const arma::vec& beta,const arma::vec& candidate,double s2,double r2,double scale,
    const arma::mat& raw_f,const arma::mat& qf,double inverse_temperature=1.0) {
  Terms t;t.raw_forward=raw_f;t.forward_precision=qf;
  arma::vec current=coordinates(alpha(j),beta);
  arma::vec next_beta=candidate.tail(beta.n_elem);
  arma::vec next_alpha=alpha;next_alpha(j)=candidate(0);
  arma::mat next_xb=xb;next_xb.col(j)=design*next_beta;
  t.raw_reverse=raw_precision(design,Y,next_alpha,pp,next_xb,j,candidate,s2,r2,inverse_temperature);
  t.reverse_precision=precision(t.raw_reverse,scale);
  arma::mat am(Y.n_rows,Y.n_cols);am.each_row()=alpha.t();
  // Existing replacement helper keeps am fixed, so place the alpha displacement
  // in the replacement predictor. All other categories remain unchanged.
  arma::vec replacement=next_xb.col(j)+(candidate(0)-alpha(j));
  if(dm_runtime::use_cached && dm_runtime::cache_replacements) {
    DMRuntimeReplacementContext context(Y,am,xb,j);
    t.likelihood_delta=context.delta(replacement);
  } else t.likelihood_delta=dm_stable_replacement_delta(Y,am,xb,j,replacement);
  t.alpha_prior_delta=-.5*(candidate(0)*candidate(0)-alpha(j)*alpha(j))/s2;
  t.beta_prior_delta=-.5*(arma::dot(next_beta,next_beta)-arma::dot(beta,beta))/r2;
  t.log_forward=dmvn_prec(candidate,current,qf);
  t.log_reverse=dmvn_prec(current,candidate,t.reverse_precision);
  t.ratio=inverse_temperature*t.likelihood_delta+t.alpha_prior_delta+t.beta_prior_delta+t.log_reverse-t.log_forward;
  return t;
}
inline Rcpp::List audit(const arma::mat& design,const arma::imat& Y,
    const arma::vec& alpha,const arma::mat& xb,int category,const arma::vec& beta,
    const arma::vec& candidate,double s2,double r2,double scale,double inverse_temperature=1.0) {
  int j=category-1;
  if(j<0||j>=static_cast<int>(Y.n_cols)||design.n_rows!=Y.n_rows||
      design.n_cols!=beta.n_elem||candidate.n_elem!=beta.n_elem+1||
      alpha.n_elem!=Y.n_cols||!candidate.is_finite()||s2<=0||r2<=0||scale<=0||
      !std::isfinite(inverse_temperature)||inverse_temperature<=0||inverse_temperature>1)
    Rcpp::stop("Invalid joint-block audit inputs");
  arma::ivec pp(Y.n_cols,arma::fill::zeros);pp(j)=beta.n_elem;
  arma::vec current=coordinates(alpha(j),beta);
  arma::mat raw=raw_precision(design,Y,alpha,pp,xb,j,current,s2,r2,inverse_temperature),qf=precision(raw,scale);
  Terms t=terms(design,Y,alpha,pp,xb,j,beta,candidate,s2,r2,scale,raw,qf,inverse_temperature);
  return Rcpp::List::create(Rcpp::Named("likelihood_delta")=t.likelihood_delta,
    Rcpp::Named("alpha_prior_delta")=t.alpha_prior_delta,Rcpp::Named("beta_prior_delta")=t.beta_prior_delta,
    Rcpp::Named("log_forward")=t.log_forward,Rcpp::Named("log_reverse")=t.log_reverse,
    Rcpp::Named("ratio")=t.ratio,Rcpp::Named("raw_forward")=t.raw_forward,
    Rcpp::Named("raw_reverse")=t.raw_reverse,Rcpp::Named("forward_precision")=t.forward_precision,
    Rcpp::Named("reverse_precision")=t.reverse_precision);
}
}
#endif
