// [[Rcpp::depends(RcppArmadillo)]]

// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include "../New_simulations/gaussian_revision_shared.h"
#include <limits>
using namespace Rcpp;

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

int sample_int(int n, int size, const arma::vec& probs, bool replace = true) {
  int pn = probs.n_elem;
  
  if (pn != n) {
    Rcpp::stop("Length of probabilities vector must match n.");
  }
  // Convert arma::vec to Rcpp::NumericVector
  Rcpp::NumericVector rcpp_probs(probs.begin(), probs.end());
  
  // Perform sampling
  return Rcpp::sample(n, size, replace, rcpp_probs)[0];
}

double RSS(const arma::vec& y_hat, const arma::vec& y) {
  
  // Compute Residual Sum of Squares (RSS)
  return arma::accu(arma::pow(y_hat - y, 2)); // arma::sum(arma::square(y - y_hat));
  
}

double marg_loglik_eff(double n, const arma::mat& XT_X, const arma::vec& XT_y, 
                       double yT_y, double a, double d, 
                       const arma::mat& V_inv,
                       double V_inv_det) {
  
  // Compute Vstar_inv = V_inv + XT_X
  arma::mat Vstar_inv = V_inv + XT_X;
  
  // Cholesky decomposition for solving and determinant calculation
  arma::mat chol_Vstar = arma::chol(Vstar_inv, "lower");
  
  // Solve for mstar using Cholesky decomposition
  arma::vec mstar = arma::solve(arma::trimatu(chol_Vstar.t()), 
                                arma::solve(arma::trimatl(chol_Vstar), XT_y));
  
  // Compute astar efficiently
  double astar = a + 0.5 * (yT_y - arma::dot(mstar, XT_y));
  
  // Update dstar
  double dstar = d + 0.5 * n;
  
  // Compute log determinant differences using Cholesky
  double logdet_Vstar = 2.0 * arma::sum(arma::trunc_log(chol_Vstar.diag()));
  
  // Compute log-marginal likelihood
  double f = -0.5 * n * (log(M_PI) + log(2)) +
    0.5 * (arma::trunc_log(V_inv_det) - logdet_Vstar) +
    d * arma::trunc_log(a) - dstar * arma::trunc_log(astar) +
    std::lgamma(dstar) - std::lgamma(d);
  
  return f;
}

// [[Rcpp::export]]
double audit_f_log_target_cpp(const arma::vec& y, const arma::mat& X,
                                     const arma::vec& gamma, double a=2.0,
                                     double d=2.0, double tau=1.0,
                                     double a_phi=1.0, double b_phi=3.0) {
  arma::uvec cols = arma::join_cols(arma::uvec({0}), arma::find(gamma > 0) + 1);
  arma::mat Xm = X.cols(cols);
  int q = cols.n_elem;
  return marg_loglik_eff(y.n_elem, Xm.t()*Xm, Xm.t()*y, arma::dot(y,y),
    a,d,arma::eye(q,q)/tau,std::pow(tau,-q)) +
    gaussrev::prior(arma::accu(gamma),gamma.n_elem,a_phi,b_phi);
}

arma::vec dbinom_cpp(const arma::vec& x, int size, double prob, bool log_density = true) {
  int n = x.n_elem;
  arma::vec densities(n); // Vector to store densities for each element of x
  
  // Check for valid input
  if (prob < 0 || prob > 1) {
    stop("Probability 'prob' must be in [0, 1].");
  }
  
  // Calculate the binomial coefficient part using lgamma for stability
  double log_binom_coeff_base = lgamma(size + 1);
  
  for (int i = 0; i < n; i++) {
    if (x(i) < 0 || x(i) > size) {
      // If x[i] is out of range, assign zero density
      densities(i) = (log_density ? -INFINITY : 0.0);
    } else {
      // Calculate log of binomial coefficient for current x
      double log_binom_coeff = log_binom_coeff_base - lgamma(x(i) + 1) - lgamma(size - x(i) + 1);
      
      // Calculate log density for current x[i]
      double log_prob_density = log_binom_coeff + x(i) * std::log(prob) + (size - x(i)) * std::log(1 - prob);
      
      // Store the result in densities vector
      densities(i) = log_density ? log_prob_density : std::exp(log_prob_density);
    }
  }
  
  return densities;
}

double F_stat(double RSS1, double RSS2, int n, int p1, int p2) {
  if (p2 == p1) return 1.0; // same model, no improvement
  double num = (RSS1 - RSS2) / (p2 - p1);
  double den = RSS2 / (n - p2);
  return num / den;
}

// Helper function to compute p-value from F-statistic
// Uses the F-distribution CDF
double f_test_pvalue(double F_stat, int df1, int df2) {
  if (F_stat <= 0) return 1.0;
  // Use boost or R::pf if available, or approximate
  // For now, using a simple approximation or assuming you have access to pf function
  // In R/Rcpp: R::pf(F_stat, df1, df2, false, false) gives upper tail probability
  // Here we assume you have a function like this available:
  return R::pf(F_stat, df1, df2, false, false); // upper tail p-value
}

arma::vec normalize_proposal(const arma::vec& log_weight, const arma::vec& gamma,
                             bool balanced_directions, double uniform_mixture) {
  const int p = gamma.n_elem;
  arma::vec prop(p, arma::fill::zeros);
  if (!balanced_directions) {
    prop = arma::exp(log_weight - log_weight.max());
    prop /= arma::sum(prop);
    if (uniform_mixture > 0.0)
      prop = (1.0 - uniform_mixture) * prop + uniform_mixture / p;
    return prop;
  }
  arma::uvec births = arma::find(gamma == 0);
  arma::uvec deaths = arma::find(gamma > 0);
  bool has_births = births.n_elem > 0, has_deaths = deaths.n_elem > 0;
  double birth_mass = has_births && has_deaths ? 0.5 : (has_births ? 1.0 : 0.0);
  double death_mass = has_births && has_deaths ? 0.5 : (has_deaths ? 1.0 : 0.0);
  if (has_births) {
    arma::vec z = arma::exp(log_weight.elem(births) - log_weight.elem(births).max());
    prop.elem(births) = birth_mass * z / arma::sum(z);
  }
  if (has_deaths) {
    arma::vec z = arma::exp(log_weight.elem(deaths) - log_weight.elem(deaths).max());
    prop.elem(deaths) = death_mass * z / arma::sum(z);
  }
  if (uniform_mixture > 0.0) {
    prop *= (1.0 - uniform_mixture);
    if (has_births) prop.elem(births) += birth_mass * uniform_mixture / births.n_elem;
    if (has_deaths) prop.elem(deaths) += death_mass * uniform_mixture / deaths.n_elem;
  }
  return prop;
}

// All one-flip F scores from one factorization of the current model.  This is
// algebraically equivalent to the legacy candidate-by-candidate regressions.
arma::vec proposal_inf_rank_one(const arma::vec& gamma,
                                const arma::vec& XT_y,
                                const arma::mat& XT_X,
                                double y_sum_sq, double RSS_null, int n,
                                double lambda, bool balanced_directions,
                                double uniform_mixture) {
  const int p = gamma.n_elem;
  const int current_k = arma::accu(gamma);
  arma::uvec selected = arma::join_cols(arma::uvec({0}), arma::find(gamma > 0) + 1);
  arma::mat gram_current = XT_X.submat(selected, selected);
  arma::mat gram_inverse;
  if (selected.n_elem >= (arma::uword)n || !arma::inv_sympd(gram_inverse, gram_current)) {
    arma::vec log_weight(p);
    for (int j = 0; j < p; ++j) {
      arma::vec candidate = gamma; candidate(j) = 1 - candidate(j);
      arma::uvec cols = arma::join_cols(arma::uvec({0}), arma::find(candidate > 0) + 1);
      double residual = gaussrev::rss(XT_X.submat(cols, cols), XT_y(cols), y_sum_sq);
      double score = gaussrev::f_score(arma::accu(candidate), n, RSS_null, residual);
      log_weight(j) = std::pow(score, lambda);
    }
    return normalize_proposal(log_weight, gamma, balanced_directions, uniform_mixture);
  }
  arma::vec beta = gram_inverse * XT_y(selected);
  double rss_current = std::max(y_sum_sq - arma::dot(XT_y(selected), beta), 1e-14);
  arma::ivec active_position(p, arma::fill::value(-1));
  for (arma::uword pos = 1; pos < selected.n_elem; ++pos)
    active_position(selected(pos) - 1) = static_cast<int>(pos);

  arma::vec log_weight(p, arma::fill::zeros);
  for (int j = 0; j < p; ++j) {
    int candidate_k = current_k + (gamma(j) > 0 ? -1 : 1);
    double rss_candidate;
    if (gamma(j) > 0) {
      int pos = active_position(j);
      rss_candidate = rss_current + beta(pos) * beta(pos) / gram_inverse(pos, pos);
    } else {
      arma::uword column = static_cast<arma::uword>(j + 1);
      arma::vec cross = XT_X.submat(selected, arma::uvec({column}));
      double denominator = XT_X(column, column) - arma::as_scalar(cross.t() * gram_inverse * cross);
      double numerator = XT_y(column) - arma::dot(cross, beta);
      rss_candidate = denominator > 1e-12
        ? rss_current - numerator * numerator / denominator : rss_current;
    }
    rss_candidate = std::max(rss_candidate, 1e-14);
    double score = gaussrev::f_score(candidate_k, n, RSS_null, rss_candidate);
    log_weight(j) = std::pow(score, lambda);
  }
  return normalize_proposal(log_weight, gamma, balanced_directions, uniform_mixture);
}

// Function to calculate proposal probabilities based on F-test p-values
// Proposal is proportional to exp((-log(pval))^lambda)
arma::vec proposal_inf(const arma::vec& gamma, const arma::vec& y, const arma::mat& X,
                       const arma::vec& XT_y, const arma::mat& XT_X,
                       double RSS_null, double lambda,
                       bool balanced_directions = false, double uniform_mixture = 0.0,
                       const std::string& proposal_engine = "legacy") {
  if (proposal_engine == "rank_one")
    return proposal_inf_rank_one(gamma, XT_y, XT_X, arma::dot(y,y), RSS_null,
                                 y.n_elem, lambda, balanced_directions, uniform_mixture);
  if (proposal_engine != "legacy") stop("proposal_engine must be legacy or rank_one");
  
  double y_sum_sq = arma::dot(y, y);              // This is y'y
  const int p = gamma.n_elem;
  const int n = y.n_elem;
  arma::vec pval(p);
  
  // Pre-compute constants
  const arma::uvec intercept_idx = {0};
  
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < p; i++) {
    // Flip bit i in gamma
    arma::vec gammanew = gamma;
    gammanew(i) = 1.0 - gammanew(i);
    
    // Select columns: intercept + active predictors
    arma::uvec active_predictors = arma::find(gammanew > 0) + 1;
    arma::uvec selected = arma::join_cols(intercept_idx, active_predictors);
    const int p_prop = selected.n_elem;
    
    // Handle edge case: if proposed model is null model (intercept only)
    if (p_prop == 1) {
      // Comparing null to null - no improvement
      pval(i) = 1.0;  // or use a large p-value
      continue;  // skip to next iteration
    }
    
    // Extract submatrices
    arma::mat XnewT_Xnew = XT_X.submat(selected, selected);
    arma::vec XnewT_y = XT_y(selected);
    
    double RSS_prop = gaussrev::rss(XnewT_Xnew, XnewT_y, y_sum_sq);
    pval(i) = std::pow(10.0, -gaussrev::f_score(p_prop - 1, n, RSS_null, RSS_prop));

  }
  
  // Weight proportional to exp((-log(pval))^lambda)
  // Store (-log(pval))^lambda, will exponentiate later
  arma::vec neg_log_pval = -arma::trunc_log(pval)/std::log(10.0);
  arma::vec log_weight = arma::pow(neg_log_pval, lambda);
  
  // Normalize using log-sum-exp trick
  // Weights are proportional to exp(log_weight)
  return normalize_proposal(log_weight, gamma, balanced_directions, uniform_mixture);
}

// Conditional probabilities for replacing one active predictor by one of a
// supplied set of inactive predictors.  The rank-one path factors the model
// with the outgoing predictor removed once, then evaluates every addition by
// an OLS rank-one identity.
arma::vec swap_probabilities_f(const arma::vec& gamma, int outgoing,
                               const arma::uvec& candidates,
                               const arma::vec& XT_y, const arma::mat& XT_X,
                               double yT_y, double RSS_null, int n,
                               double lambda, const std::string& engine,
                               const std::string& weight_mode) {
  if (candidates.n_elem == 0) return arma::vec();
  if (weight_mode == "uniform")
    return arma::vec(candidates.n_elem, arma::fill::value(1.0 / candidates.n_elem));
  if (weight_mode != "informed") stop("move_weight_mode must be informed or uniform");
  arma::vec log_weight(candidates.n_elem, arma::fill::zeros);
  arma::vec base_gamma = gamma;
  base_gamma(outgoing) = 0;
  arma::uvec base = join_cols(arma::uvec({0}), arma::find(base_gamma > 0) + 1);
  arma::mat base_inverse;
  arma::vec base_beta;
  double rss_base = 0.0;
  bool rank_one_ok = false;
  if (engine == "rank_one") {
    arma::mat gram = XT_X.submat(base, base);
    rank_one_ok = base.n_elem < (arma::uword)n && arma::inv_sympd(base_inverse, gram);
    if (!rank_one_ok) base_inverse = arma::pinv(gram);
    base_beta = base_inverse * XT_y(base);
    rss_base = std::max(yT_y - arma::dot(XT_y(base), base_beta), 1e-14);
  } else if (engine != "legacy") stop("proposal_engine must be legacy or rank_one");
  const int candidate_k = arma::accu(gamma);
  const int p_prop = candidate_k + 1;
  const int df2 = std::max(n - p_prop, 1);
  for (arma::uword j = 0; j < candidates.n_elem; ++j) {
    double rss;
    if (engine == "rank_one" && rank_one_ok) {
      arma::uword column = candidates(j) + 1;
      arma::vec cross = XT_X.submat(base, arma::uvec({column}));
      double denominator = XT_X(column,column) - arma::as_scalar(cross.t()*base_inverse*cross);
      double numerator = XT_y(column) - arma::dot(cross,base_beta);
      rss = denominator > 1e-12 ? rss_base - numerator*numerator/denominator : rss_base;
    } else {
      arma::vec candidate_gamma = base_gamma;
      candidate_gamma(candidates(j)) = 1;
      arma::uvec selected = join_cols(arma::uvec({0}), arma::find(candidate_gamma > 0) + 1);
      arma::mat gram = XT_X.submat(selected,selected);
      arma::vec cross_y = XT_y(selected);
      rss = gaussrev::rss(gram, cross_y, yT_y);
    }
    rss = std::max(rss,1e-14);
    double score = gaussrev::f_score(candidate_k, n, RSS_null, rss);
    log_weight(j) = std::pow(score,lambda);
  }
  arma::vec weights = arma::exp(log_weight-log_weight.max());
  return weights/arma::sum(weights);
}

// [[Rcpp::export]]
arma::vec audit_f_swap_proposal_cpp(const arma::vec& y, const arma::mat& X,
                                    const arma::vec& gamma, int outgoing_one_based,
                                    const arma::uvec& candidates_one_based,
                                    double lambda=1.0,
                                    std::string proposal_engine="legacy",
                                    std::string move_weight_mode="informed") {
  return swap_probabilities_f(gamma,outgoing_one_based-1,candidates_one_based-1,
    X.t()*y,X.t()*X,arma::dot(y,y),arma::accu(arma::square(y-arma::mean(y))),
    y.n_elem,lambda,proposal_engine,move_weight_mode);
}

// [[Rcpp::export]]
arma::vec audit_f_proposal_cpp(const arma::vec& y, const arma::mat& X,
                               const arma::vec& gamma, double lambda = 0.7,
                               bool balanced_directions = false,
                               double uniform_mixture = 0.0,
                               std::string proposal_engine = "legacy") {
  arma::mat XT_X = X.t() * X;
  arma::vec XT_y = X.t() * y;
  double RSS_null = arma::accu(arma::square(y - arma::mean(y)));
  return proposal_inf(gamma, y, X, XT_y, XT_X, RSS_null, lambda,
                      balanced_directions, uniform_mixture, proposal_engine);
}

// [[Rcpp::export]]
List var_sel_distance_cpp(arma::vec y, arma::mat X, std::vector<arma::vec> neigh, 
                          arma::vec gamma, double phi,
                          double a, double d, double tau, double a_phi, double b_phi,
                          int runs, int burn, int thin = 1, int start_adapt = 100, int end_adapt = 1000, 
                          double lambda_bd = 1.0, double lambda_move = 1.0, 
                          double trgt_acpt = 0.234, bool adapt = false,
                          bool move_step = false,
                          std::string adaptation_mode = "legacy",
                          bool balanced_directions = false,
                          double uniform_mixture = 0.0,
                          double target_effective_fraction = 0.025,
                          double adaptation_gain = 0.5,
                          double adaptation_decay = 0.6,
                          double adaptation_offset = 100.0,
                          std::string proposal_engine = "legacy",
                          std::string move_weight_mode = "informed") {
  
  int n = y.n_elem;
  int p = X.n_cols - 1;
  if (adaptation_mode != "legacy" && adaptation_mode != "diversity")
    stop("adaptation_mode must be legacy or diversity");
  if (uniform_mixture < 0.0 || uniform_mixture >= 1.0)
    stop("uniform_mixture must be in [0,1)");
  if (adaptation_decay <= 0.5 || adaptation_decay > 1.0)
    stop("adaptation_decay must be in (0.5,1]");
  if (proposal_engine != "legacy" && proposal_engine != "rank_one")
    stop("proposal_engine must be legacy or rank_one");
  if (move_weight_mode != "informed" && move_weight_mode != "uniform")
    stop("move_weight_mode must be informed or uniform");
  bool diversity_adaptation = adaptation_mode == "diversity";
  if (diversity_adaptation) {
    balanced_directions = true;
    if (uniform_mixture == 0.0) uniform_mixture = 0.02;
  }
  double target_effective_neighbors = std::max(5.0, target_effective_fraction * p);
  
  double RSS_null = arma::accu(arma::pow(arma::mean(y) - y, 2)); 
  
  // arma::vec gamma(p, arma::fill::ones);
  arma::mat XT_X = X.t() * X;
  arma::vec XT_y = X.t() * y;
  double yT_y = arma::as_scalar(y.t() * y);
  
  arma::uvec selected = join_cols(arma::uvec({0}), (arma::find(gamma > 0) + 1));
  
  int qinit = arma::sum(gamma);
  
  arma::mat Xinit = X.cols(selected);
  arma::mat XinitT_Xinit = XT_X.submat(selected, selected);
  arma::vec XinitT_y = XT_y(selected); 
  
  // Initial log likelihood and prior
  double loglik = marg_loglik_eff(n, XinitT_Xinit, XinitT_y, yT_y, a, d, arma::eye(qinit + 1, qinit + 1) / tau, std::pow(tau, -(qinit + 1)));
  double logprior = gaussrev::prior(arma::accu(gamma), gamma.n_elem, a_phi, b_phi);
  
  int MCMC;
  if(thin == 0){
    MCMC = floor(runs-burn);
  }else{
    MCMC = floor((runs-burn)/thin);
  }
  
  
  // Tracking variables  
  int acc_bd_window = 0;               // Acceptances in current window
  int iter_in_window = 0;              // Iterations in current window
  double prev_window_rate = 0.5;       // Previous window's acceptance rate
  
  // Add these variables at function start:
  double log_lambda_bd = 0.0;          // Work in log-space: log(λ_BD)
  double prev_accept_rate = 0.0;       // α^(j-2)
  double curr_accept_rate = 0.0;       // α^(j-1)
  double prev_log_lambda = 0.0;        // log(λ^(j-2))
  
  int update_counter = 0;              // Counts updates (j in the formula)
  double decay_power = 0.75;            // Robbins-Monro decay exponent
  
  // Tracking for current iteration
  int acc_bd_iter = 0;                 // Acceptances this iteration window
  int iter_counter = 0;                // Iterations since last update
  double proposal_ess_window = 0.0;
  
  int window_size = 25;                // How often to update
  int adaptation_start = start_adapt;
  int adaptation_end = end_adapt;
  
  // Bounds
  double log_lambda_min = std::log(0.05);
  double log_lambda_max = std::log(10.0);
  
  // Initialize
  log_lambda_bd = std::log(lambda_bd);
  
  arma::mat GAMMAVEC(MCMC, p, arma::fill::zeros);
  arma::vec phivec(MCMC, arma::fill::zeros);
  arma::vec indvec(MCMC, arma::fill::zeros);
  arma::vec loglikvec(MCMC, arma::fill::zeros);
  arma::vec accratevec(runs, arma::fill::zeros);
  arma::vec lambdavec(runs, arma::fill::zeros);
  arma::vec proposal_ess_vec(runs, arma::fill::zeros);
  
  int acc_bd = 0;
  int acc_move = 0;
  int move_attempts = 0;
  
  arma::vec prop = proposal_inf(gamma, y, X, XT_y, XT_X, RSS_null, lambda_bd,
                                balanced_directions, uniform_mixture, proposal_engine);
  
  Rcout << "### Informed proposal: F statistics" << std::endl;
  Rcout << "### lambda_bd = " << lambda_bd << std::endl;
  Rcout << "### lambda_move = " << lambda_move << std::endl;
  Rcout << " " << std::endl;
  
  double acc_rate_bd = 0.0;
  double acc_rate_move = 0.0;
  
  int n_move_attempts = 1;
  
  for (int iter = 0; iter < runs; iter++) {
    double proposal_ess = 1.0 / arma::accu(arma::square(prop));
    
    double progress = static_cast<double>(iter * 100) / runs;
    acc_rate_bd = static_cast<double>(acc_bd) / static_cast<double>(iter + 1);
    acc_rate_move = static_cast<double>(acc_move) / static_cast<double>(iter + 1);
    
    if (progress >= 10.0 && progress <= 100.0 && fmod(progress, 10.0) == 0.0) {
      Rcout << "### Progress: " << progress << " %" << std::endl;
      Rcout << "### lambda_bd = " << lambda_bd << std::endl;
      Rcout << "### lambda_move = " << lambda_move << std::endl;
      Rcout << "### Acceptance rate birth and death step: " << acc_rate_bd << std::endl;
      Rcout << "### Acceptance rate move step: " << acc_rate_move << std::endl;
      Rcout << " " << std::endl;
    }
    
    ////////////////////////
    // Birth and death step//
    ////////////////////////
    
    int ind = sample_int(p, 1, prop) - 1;
    
    arma::vec gammanew = gamma;
    gammanew(ind) = std::abs(gamma(ind) - 1);
    arma::uvec selected_new = join_cols(arma::uvec({0}), (arma::find(gammanew > 0) + 1));
    
    int qnew = arma::sum(gammanew);
    double logprior_new = gaussrev::prior(arma::accu(gammanew), gammanew.n_elem, a_phi, b_phi);
    
    arma::mat Xnew = X.cols(selected_new);
    arma::mat XnewT_Xnew = XT_X.submat(selected_new, selected_new);
    arma::vec XnewT_y = XT_y(selected_new);
    
    double loglik_new = marg_loglik_eff(n, XnewT_Xnew, XnewT_y, yT_y, a, d, 
                                        arma::eye(qnew + 1, qnew + 1) / tau, 
                                        std::pow(tau, -(qnew + 1)));
    double u = arma::randu();
    
    arma::vec prop_new = proposal_inf(gammanew, y, X, XT_y, XT_X, RSS_null, lambda_bd,
                                      balanced_directions, uniform_mixture, proposal_engine);
    
    double log_ratio = (loglik_new - loglik) + (logprior_new - logprior) + 
      (arma::trunc_log(prop_new(ind)) - arma::trunc_log(prop(ind)));
    double log_acpt = std::min(0.0, log_ratio);
    
    bool accepted = (log(u) < log_acpt);
    
    if (accepted) {
      gamma = gammanew;
      selected = selected_new;
      loglik = loglik_new;
      logprior = logprior_new;
      prop = prop_new;
      acc_bd += 1;
      acc_bd_iter += 1;
    }
    
    iter_counter += 1;
    proposal_ess_window += proposal_ess;
    
    ////////////////////////
    // ROBBINS-MONRO ADAPTATION (matching reference formula)
    ////////////////////////
    
    if (iter >= adaptation_start && (diversity_adaptation || iter < adaptation_end) &&
        iter_counter >= window_size && adapt == true) {
      
      // Compute current window acceptance rate
      curr_accept_rate = static_cast<double>(acc_bd_iter) / window_size;
      
      if (diversity_adaptation) {
        double mean_effective = proposal_ess_window / iter_counter;
        double step_size = adaptation_gain *
          std::pow(adaptation_offset + update_counter, -adaptation_decay);
        log_lambda_bd += step_size *
          (std::log(mean_effective) - std::log(target_effective_neighbors));
        log_lambda_bd = std::max(std::log(0.01), std::min(log_lambda_bd, std::log(2.0)));
        lambda_bd = std::exp(log_lambda_bd);
        prop = proposal_inf(gamma, y, X, XT_y, XT_X, RSS_null, lambda_bd,
                            balanced_directions, uniform_mixture, proposal_engine);
      } else if (update_counter > 0) {
        
        // Step size: decays as (j-1)^{-decay_power}
        double step_size = std::pow(static_cast<double>(update_counter), -decay_power);
        
        // Acceptance change: Acc^(j-1) - Acc^(j-2)
        double acc_change = curr_accept_rate - prev_accept_rate;
        
        // Direction of last lambda change: I[λ^(j-1) > λ^(j-2)] - I[λ^(j-1) ≤ λ^(j-2)]
        double lambda_direction = (log_lambda_bd > prev_log_lambda) ? 1.0 : -1.0;
        
        // Gradient signal: positive if both increased together or both decreased together
        double gradient = acc_change * lambda_direction;
        
        // Robbins-Monro update in log-space
        double log_lambda_new = log_lambda_bd + step_size * gradient;
        
        // Store previous for next iteration
        prev_log_lambda = log_lambda_bd;
        
        // Update with bounds
        log_lambda_bd = std::max(log_lambda_min, std::min(log_lambda_new, log_lambda_max));
        
        // Convert back to regular space
        lambda_bd = std::exp(log_lambda_bd);

        // Refresh the cached proposal after changing lambda; otherwise the
        // next transition uses weights from the previous adaptation value.
        prop = proposal_inf(gamma, y, X, XT_y, XT_X, RSS_null, lambda_bd,
                            balanced_directions, uniform_mixture, proposal_engine);
        
      } else {
        // First update - just store current state
        prev_log_lambda = log_lambda_bd;
      }
      
      // Update counters
      prev_accept_rate = curr_accept_rate;
      update_counter += 1;
      acc_bd_iter = 0;
      iter_counter = 0;
      proposal_ess_window = 0.0;
      
    }
    
    accratevec(iter) = static_cast<double>(acc_bd) / static_cast<double>(iter + 1);
    lambdavec(iter) = lambda_bd;
    proposal_ess_vec(iter) = proposal_ess;
    
    /////////////////
    // Move step //
    /////////////////
    
    if (move_step == true){
      bool move_accepted_this_iter = false;  // Track if ANY move accepted this iteration
    
    if (selected.n_elem > 1) {
      arma::uvec active_candidates = selected.rows(1, selected.n_elem - 1) - 1;
      if (active_candidates.n_elem > 0) {
        for (int attempt = 0; attempt < n_move_attempts; ++attempt) {
          
          // [Your existing move step code...]
          arma::uvec has_inactive_idx;
          has_inactive_idx.reset();
          for (size_t k = 0; k < active_candidates.n_elem; ++k) {
            int a = active_candidates(k);
            arma::vec neigh_a = neigh[a];
            bool ok=false;
            for (size_t t=0;t<neigh_a.n_elem;++t) if (gamma[neigh_a(t)-1]==0) { ok=true; break; }
            if (ok) has_inactive_idx.insert_rows(has_inactive_idx.n_elem, arma::uvec({(arma::uword)k}));
          }
          if (has_inactive_idx.n_elem == 0) break;
          
          int kselpos = arma::randi(arma::distr_param(0, has_inactive_idx.n_elem-1));
          int move_idx = active_candidates(has_inactive_idx(kselpos));
          
          arma::vec neigh_a = neigh[move_idx];
          std::vector<int> inactive_neigh;
          for (size_t t=0;t<neigh_a.n_elem;++t) {
            int cand = (int)neigh_a(t)-1;
            if (gamma(cand)==0) inactive_neigh.push_back(cand);
          }
          if (inactive_neigh.size()==0) continue;
          move_attempts += 1;
          
          // Score inactive neighbors using lambda_move
          int m = inactive_neigh.size();
          arma::uvec inactive_u(m);
          for (int j=0;j<m;++j) inactive_u(j)=inactive_neigh[j];
          arma::vec w = swap_probabilities_f(gamma,move_idx,inactive_u,XT_y,XT_X,
            yT_y,RSS_null,n,lambda_move,proposal_engine,move_weight_mode);
          
          double ru = arma::randu();
          double csum=0; int jpick=0;
          for (int j=0;j<m;++j){ csum += w(j); if (ru <= csum) { jpick=j; break; } }
          int rand_neigh_idx = inactive_neigh[jpick];
          
          arma::vec gammanew = gamma;
          gammanew(move_idx)=0;
          gammanew(rand_neigh_idx)=1;
          arma::uvec selected_new = join_cols(arma::uvec({0}), (arma::find(gammanew > 0) + 1));
          int qnew = arma::sum(gammanew);
          arma::mat XnewT_Xnew = XT_X.submat(selected_new, selected_new);
          arma::vec XnewT_y = XT_y(selected_new);
          double logprior_new = gaussrev::prior(arma::accu(gammanew), gammanew.n_elem, a_phi, b_phi);
          double loglik_new = marg_loglik_eff(n, XnewT_Xnew, XnewT_y, yT_y, a, d,
                                              arma::eye(qnew+1,qnew+1)/tau,
                                              std::pow(tau,-(qnew+1)));
          
          double pf = (1.0 / (double)has_inactive_idx.n_elem) * w(jpick);
          
          // [Reverse probability calculation - your existing code...]
          arma::uvec active_rev = selected_new.rows(1, selected_new.n_elem - 1) - 1;
          std::vector<int> act_rev_with_inactive;
          for (size_t rr=0; rr<active_rev.n_elem; ++rr) {
            int a = active_rev(rr);
            arma::vec neigh_a2 = neigh[a];
            bool ok=false;
            for (size_t t=0;t<neigh_a2.n_elem;++t){
              int cand = (int)neigh_a2(t)-1;
              if (gammanew[cand]==0) { ok=true; break; }
            }
            if (ok) act_rev_with_inactive.push_back(a);
          }
          
          bool reverse_possible = false;
          for (size_t i=0; i<act_rev_with_inactive.size(); ++i) {
            if (act_rev_with_inactive[i] == rand_neigh_idx) {
              reverse_possible = true;
              break;
            }
          }
          if (!reverse_possible) continue;
          
          arma::vec neigh_rev_vec = neigh[rand_neigh_idx];
          std::vector<int> inactive_rev;
          for (size_t t=0;t<neigh_rev_vec.n_elem;++t){
            int cand = (int)neigh_rev_vec(t)-1;
            if (gammanew[cand]==0) inactive_rev.push_back(cand);
          }
          if (inactive_rev.size() == 0) continue;
          
          bool move_idx_in_reverse = false;
          int jrev_index = -1;
          for (size_t j=0; j<inactive_rev.size(); ++j) {
            if (inactive_rev[j] == move_idx) {
              move_idx_in_reverse = true;
              jrev_index = j;
              break;
            }
          }
          if (!move_idx_in_reverse) continue;
          
          int mr = inactive_rev.size();
          arma::uvec inactive_rev_u(mr);
          for (int j=0;j<mr;++j) inactive_rev_u(j)=inactive_rev[j];
          arma::vec wr = swap_probabilities_f(gammanew,rand_neigh_idx,inactive_rev_u,
            XT_y,XT_X,yT_y,RSS_null,n,lambda_move,proposal_engine,move_weight_mode);
          
          double pr = (1.0 / (double)act_rev_with_inactive.size()) * wr(jrev_index);
          
          double proposal_log_ratio = std::log(pr) - std::log(pf);
          double log_ratio = (loglik_new - loglik) + (logprior_new - logprior) + proposal_log_ratio;
          
          // Track this attempt
          // move_attempts_window += 1;
          
          if (std::log(arma::randu()) < std::min(0.0, log_ratio)) {
            // ACCEPT
            gamma = gammanew;
            selected = selected_new;
            loglik = loglik_new;
            logprior = logprior_new;
            prop = proposal_inf(gamma, y, X, XT_y, XT_X, RSS_null, lambda_bd,
                                balanced_directions, uniform_mixture, proposal_engine);
            acc_move += 1;
            // acc_move_window += 1;
            move_accepted_this_iter = true;
            break; // Stop further attempts this iteration
          }
        } // end attempts loop
      }
    }
    }
    
    ////////////////
    // Auxiliary conditional phi draw; the MH target integrates phi out. //
    ////////////////
    double sum_gamma = arma::accu(gamma);
    phi = R::rbeta(a_phi + sum_gamma, b_phi + p - sum_gamma);
    
    if((fmod(iter, thin) == 0.0) && (iter >= burn)){
      GAMMAVEC.row(iter - burn) = gamma.t();
      phivec(iter - burn) = phi;
      indvec(iter-burn) = ind;
      loglikvec(iter-burn) = loglik + logprior;
    }
    
  }
  
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("gamma") = GAMMAVEC,
    Rcpp::Named("prop") = indvec,
    Rcpp::Named("phi") = phivec,
    Rcpp::Named("loglik") = loglikvec,
    Rcpp::Named("acc_rate_bd") = acc_rate_bd,
    Rcpp::Named("acc_rate_move") = acc_rate_move,
    Rcpp::Named("move_attempts") = move_attempts,
    Rcpp::Named("move_acceptances") = acc_move,
    Rcpp::Named("conditional_move_acceptance") = move_attempts > 0 ?
      static_cast<double>(acc_move) / move_attempts : 0.0,
    Rcpp::Named("acc_rate_bd_chain") = accratevec,
    Rcpp::Named("lambda_chain") = lambdavec,
    Rcpp::Named("proposal_effective_neighbors") = proposal_ess_vec,
    Rcpp::Named("adaptation_mode") = adaptation_mode,
    Rcpp::Named("proposal_engine") = proposal_engine,
    Rcpp::Named("move_weight_mode") = move_weight_mode,
    Rcpp::Named("target_effective_neighbors") = target_effective_neighbors
  );
  output["final_proposal"] = prop;
  output["final_gamma"] = gamma;
  output["target_prior"] = "collapsed_beta_binomial";
  output["balanced_directions"] = balanced_directions;
  output["uniform_mixture"] = uniform_mixture;
  return output;
  
}
