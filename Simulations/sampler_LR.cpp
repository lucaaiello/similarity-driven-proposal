// [[Rcpp::depends(RcppArmadillo)]]

// #include <Rcpp.h>
#include <RcppArmadillo.h>
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

double LR_test_pvalue(double lik_ratio_stat, int df) {
  if (lik_ratio_stat <= 0) return 1.0;
  return R::pchisq(lik_ratio_stat, df, false, false); // upper tail p-value
}

arma::vec proposal_inf(const arma::vec& gamma, const arma::vec& y, const arma::mat& X,
                       const arma::vec& XT_y, const arma::mat& XT_X,
                       double RSS_null, double lambda) {
  
  double y_sum_sq = arma::dot(y, y);              // This is y'y
  
  const int p = gamma.n_elem;
  const int n = y.n_elem;
  arma::vec pval(p);
  
  const arma::uvec intercept_idx = {0};
  
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < p; ++i) {
    // Flip bit i in gamma
    arma::vec gammanew = gamma;
    gammanew(i) = 1.0 - gammanew(i);
    
    // Select columns: intercept + active predictors
    arma::uvec active_predictors = arma::find(gammanew > 0) + 1;
    arma::uvec selected = arma::join_cols(intercept_idx, active_predictors);
    const int p_prop = selected.n_elem;
    
    // Extract submatrices
    arma::mat XnewT_Xnew = XT_X.submat(selected, selected);
    arma::vec XnewT_y = XT_y(selected);
    
    // Solve normal equations
    arma::vec betahat = arma::solve(XnewT_Xnew, XnewT_y);
    
    // Compute RSS for proposed model
    double RSS_prop = y_sum_sq - arma::dot(XnewT_y, betahat);  // This is y'y - β'(X'y)
    // RSS_prop = std::max(RSS_prop, 1e-12);
    
    // LRT statistic
    double lik_ratio_stat = n * (std::log(RSS_null) - std::log(RSS_prop));
    // int df = p_prop - 1;               // number of covariates
    // df = std::max(df, 1);
    // 
    // // upper-tail p-value, floor very small values
    // pval(i) = LR_test_pvalue(lik_ratio_stat, df);
    
    // // Per-parameter LRT
    // double lik_ratio_per_param = lik_ratio_stat / (p_prop - 1);
    // 
    // // Use df=1 to interpret as "average parameter quality"
    // int df = 1;
    // pval(i) = LR_test_pvalue(lik_ratio_per_param, df);
    
    // OPTION 3: Power-adjusted LRT
    double alpha = 0.5;
    double lik_ratio_adjusted = lik_ratio_stat / std::pow(p_prop - 1, alpha);
    int df = 1;
    pval(i) = LR_test_pvalue(lik_ratio_adjusted, df);
    
    // OPTION 4: BIC-style penalty
    // double lik_ratio_penalized = lik_ratio_stat - std::log(n) * (p_prop - 1);
    // lik_ratio_penalized = std::max(lik_ratio_penalized, 0.0);
    // int df = 1;
    // pval(i) = LR_test_pvalue(lik_ratio_penalized, df);
    
  }
  
  arma::vec neg_log_pval = -arma::trunc_log(pval)/std::log(10.0);
  arma::vec log_weight = arma::pow(neg_log_pval, lambda);
  
  // Normalize using log-sum-exp trick
  // Weights are proportional to exp(log_weight)
  double max_weight = log_weight.max();
  arma::vec prop = arma::exp(log_weight - max_weight);
  prop /= arma::sum(prop);
  
  return prop;
}


// [[Rcpp::export]]
List var_sel_distance_cpp(arma::vec y, arma::mat X, std::vector<arma::vec> neigh, 
                          arma::vec gamma, double phi,
                          double a, double d, double tau, double a_phi, double b_phi,
                          int runs, int burn, int thin = 1, int start_adapt = 100, int end_adapt = 1000, 
                          double lambda_bd = 1.0, double lambda_move = 1.0, 
                          double trgt_acpt = 0.234, bool adapt = false,
                          bool move_step = true) {
  
  int n = y.n_elem;
  int p = X.n_cols - 1;
  
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
  double logprior = arma::sum(dbinom_cpp(gamma, 1, phi));
  
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
  
  int acc_bd = 0;
  int acc_move = 0;
  
  arma::vec prop = proposal_inf(gamma, y, X, XT_y, XT_X, RSS_null, lambda_bd);
  
  Rcout << "### Informed proposal: F statistics" << std::endl;
  Rcout << "### lambda_bd = " << lambda_bd << std::endl;
  Rcout << "### lambda_move = " << lambda_move << std::endl;
  Rcout << " " << std::endl;
  
  double acc_rate_bd = 0.0;
  double acc_rate_move = 0.0;
  
  int n_move_attempts = 1;
  
  for (int iter = 0; iter < runs; iter++) {
    
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
    double logprior_new = arma::sum(dbinom_cpp(gammanew, 1, phi));
    
    arma::mat Xnew = X.cols(selected_new);
    arma::mat XnewT_Xnew = XT_X.submat(selected_new, selected_new);
    arma::vec XnewT_y = XT_y(selected_new);
    
    double loglik_new = marg_loglik_eff(n, XnewT_Xnew, XnewT_y, yT_y, a, d, 
                                        arma::eye(qnew + 1, qnew + 1) / tau, 
                                        std::pow(tau, -(qnew + 1)));
    double u = arma::randu();
    
    arma::vec prop_new = proposal_inf(gammanew, y, X, XT_y, XT_X, RSS_null, lambda_bd);
    
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
    
    ////////////////////////
    // ROBBINS-MONRO ADAPTATION (matching reference formula)
    ////////////////////////
    
    if (iter >= adaptation_start && iter < adaptation_end &&
        iter_counter >= window_size && adapt == true) {
      
      // Compute current window acceptance rate
      curr_accept_rate = static_cast<double>(acc_bd_iter) / window_size;
      
      // Only update if we have a previous rate (after first window)
      if (update_counter > 0) {
        
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
        
      } else {
        // First update - just store current state
        prev_log_lambda = log_lambda_bd;
      }
      
      // Update counters
      prev_accept_rate = curr_accept_rate;
      update_counter += 1;
      acc_bd_iter = 0;
      iter_counter = 0;
      
    }
    
    accratevec(iter) = static_cast<double>(acc_bd) / static_cast<double>(iter + 1);
    lambdavec(iter) = lambda_bd;
    
    /////////////////
    // Move step //
    /////////////////
    
    if (move_step == true){
      bool move_accepted_this_iter = false;  // Track if ANY move accepted this iteration
    
    if (selected.n_elem > 1) {
      arma::uvec active_candidates = selected.rows(1, selected.n_elem - 1) - 1;
      if (active_candidates.n_elem > 0) {
        for (int attempt = 0; attempt < n_move_attempts; ++attempt) {
          
          // Find active variables with inactive neighbors
          arma::uvec has_inactive_idx;
          has_inactive_idx.reset();
          for (size_t k = 0; k < active_candidates.n_elem; ++k) {
            int a = active_candidates(k);
            arma::vec neigh_a = neigh[a];
            bool ok = false;
            for (size_t t = 0; t < neigh_a.n_elem; ++t) {
              if (gamma[neigh_a(t) - 1] == 0) { 
                ok = true; 
                break; 
              }
            }
            if (ok) has_inactive_idx.insert_rows(has_inactive_idx.n_elem, arma::uvec({(arma::uword)k}));
          }
          if (has_inactive_idx.n_elem == 0) break;
          
          // Uniformly select an active variable with inactive neighbors
          int kselpos = arma::randi(arma::distr_param(0, has_inactive_idx.n_elem - 1));
          int move_idx = active_candidates(has_inactive_idx(kselpos));
          
          // Get inactive neighbors
          arma::vec neigh_a = neigh[move_idx];
          std::vector<int> inactive_neigh;
          for (size_t t = 0; t < neigh_a.n_elem; ++t) {
            int cand = (int)neigh_a(t) - 1;
            if (gamma(cand) == 0) inactive_neigh.push_back(cand);
          }
          if (inactive_neigh.size() == 0) continue;
          
          // Score inactive neighbors using likelihood ratio test
          int m = inactive_neigh.size();
          arma::vec log_weight(m);
          
          for (int j = 0; j < m; ++j) {
            int cand0 = inactive_neigh[j];
            arma::vec gamm_tmp = gamma;
            gamm_tmp(move_idx) = 0;
            gamm_tmp(cand0) = 1;
            arma::uvec sel_tmp = join_cols(arma::uvec({0}), (arma::find(gamm_tmp > 0) + 1));
            int p_prop = sel_tmp.n_elem;
            
            arma::mat XtX_tmp = XT_X.submat(sel_tmp, sel_tmp);
            arma::vec Xty_tmp = XT_y(sel_tmp);
            arma::vec beta_tmp = arma::solve(XtX_tmp, Xty_tmp);
            double RSS_tmp = yT_y - arma::dot(Xty_tmp, beta_tmp);
            RSS_tmp = std::max(RSS_tmp, 1e-12);
            
            // LRT statistic
            double lik_ratio_stat = n * (std::log(RSS_null) - std::log(RSS_tmp));
            
            // Power-adjusted LRT (matching proposal_inf)
            double alpha = 0.5;
            double lik_ratio_adjusted = lik_ratio_stat / std::pow(p_prop - 1, alpha);
            int df = 1;
            double pval = LR_test_pvalue(lik_ratio_adjusted, df);
            pval = std::max(pval, 1e-300);
            
            double neg_log_pval = -std::log10(pval);
            log_weight(j) = std::pow(neg_log_pval, lambda_move);
          }
          
          // Normalize using log-sum-exp trick
          double max_weight = log_weight.max();
          arma::vec w = arma::exp(log_weight - max_weight);
          w /= arma::sum(w);
          
          // Sample inactive neighbor
          double ru = arma::randu();
          double csum = 0; 
          int jpick = 0;
          for (int j = 0; j < m; ++j) { 
            csum += w(j); 
            if (ru <= csum) { 
              jpick = j; 
              break; 
            } 
          }
          int rand_neigh_idx = inactive_neigh[jpick];
          
          // Construct proposed model
          arma::vec gammanew = gamma;
          gammanew(move_idx) = 0;
          gammanew(rand_neigh_idx) = 1;
          arma::uvec selected_new = join_cols(arma::uvec({0}), (arma::find(gammanew > 0) + 1));
          int qnew = arma::sum(gammanew);
          arma::mat XnewT_Xnew = XT_X.submat(selected_new, selected_new);
          arma::vec XnewT_y = XT_y(selected_new);
          double logprior_new = arma::sum(dbinom_cpp(gammanew, 1, phi));
          double loglik_new = marg_loglik_eff(n, XnewT_Xnew, XnewT_y, yT_y, a, d,
                                              arma::eye(qnew + 1, qnew + 1) / tau,
                                              std::pow(tau, -(qnew + 1)));
          
          double pf = (1.0 / (double)has_inactive_idx.n_elem) * w(jpick);
          
          // Compute reverse probability
          arma::uvec active_rev = selected_new.rows(1, selected_new.n_elem - 1) - 1;
          std::vector<int> act_rev_with_inactive;
          for (size_t rr = 0; rr < active_rev.n_elem; ++rr) {
            int a = active_rev(rr);
            arma::vec neigh_a2 = neigh[a];
            bool ok = false;
            for (size_t t = 0; t < neigh_a2.n_elem; ++t) {
              int cand = (int)neigh_a2(t) - 1;
              if (gammanew[cand] == 0) { 
                ok = true; 
                break; 
              }
            }
            if (ok) act_rev_with_inactive.push_back(a);
          }
          
          bool reverse_possible = false;
          for (size_t i = 0; i < act_rev_with_inactive.size(); ++i) {
            if (act_rev_with_inactive[i] == rand_neigh_idx) {
              reverse_possible = true;
              break;
            }
          }
          if (!reverse_possible) continue;
          
          // Find inactive neighbors in reverse direction
          arma::vec neigh_rev_vec = neigh[rand_neigh_idx];
          std::vector<int> inactive_rev;
          for (size_t t = 0; t < neigh_rev_vec.n_elem; ++t) {
            int cand = (int)neigh_rev_vec(t) - 1;
            if (gammanew[cand] == 0) inactive_rev.push_back(cand);
          }
          if (inactive_rev.size() == 0) continue;
          
          bool move_idx_in_reverse = false;
          int jrev_index = -1;
          for (size_t j = 0; j < inactive_rev.size(); ++j) {
            if (inactive_rev[j] == move_idx) {
              move_idx_in_reverse = true;
              jrev_index = j;
              break;
            }
          }
          if (!move_idx_in_reverse) continue;
          
          // Score reverse neighbors using likelihood ratio test
          int mr = inactive_rev.size();
          arma::vec log_weight_rev(mr);
          
          for (int j = 0; j < mr; ++j) {
            int cand0 = inactive_rev[j];
            arma::vec gamm_tmp = gammanew;
            gamm_tmp(rand_neigh_idx) = 0;
            gamm_tmp(cand0) = 1;
            arma::uvec sel_tmp = join_cols(arma::uvec({0}), (arma::find(gamm_tmp > 0) + 1));
            int p_prop = sel_tmp.n_elem;
            
            arma::mat XtX_tmp = XT_X.submat(sel_tmp, sel_tmp);
            arma::vec Xty_tmp = XT_y(sel_tmp);
            arma::vec beta_tmp = arma::solve(XtX_tmp, Xty_tmp);
            double RSS_tmp = yT_y - arma::dot(Xty_tmp, beta_tmp);
            RSS_tmp = std::max(RSS_tmp, 1e-12);
            
            // LRT statistic
            double lik_ratio_stat = n * (std::log(RSS_null) - std::log(RSS_tmp));
            
            // Power-adjusted LRT (matching proposal_inf)
            double alpha = 0.5;
            double lik_ratio_adjusted = lik_ratio_stat / std::pow(p_prop - 1, alpha);
            int df = 1;
            double pval = LR_test_pvalue(lik_ratio_adjusted, df);
            pval = std::max(pval, 1e-300);
            
            double neg_log_pval = -std::log10(pval);
            log_weight_rev(j) = std::pow(neg_log_pval, lambda_move);
          }
          
          // Normalize reverse weights
          double max_weight_rev = log_weight_rev.max();
          arma::vec wr = arma::exp(log_weight_rev - max_weight_rev);
          wr /= arma::sum(wr);
          
          double pr = (1.0 / (double)act_rev_with_inactive.size()) * wr(jrev_index);
          
          // Metropolis-Hastings ratio
          double proposal_log_ratio = std::log(pr) - std::log(pf);
          double log_ratio = (loglik_new - loglik) + (logprior_new - logprior) + proposal_log_ratio;
          
          if (std::log(arma::randu()) < std::min(0.0, log_ratio)) {
            // ACCEPT
            gamma = gammanew;
            selected = selected_new;
            loglik = loglik_new;
            logprior = logprior_new;
            prop = proposal_inf(gamma, y, X, XT_y, XT_X, RSS_null, lambda_bd);
            acc_move += 1;
            move_accepted_this_iter = true;
            break; // Stop further attempts this iteration
          }
        } // end attempts loop
      }
    }
    }
    
    ////////////////
    // Update phi //
    ////////////////
    double sum_gamma = arma::accu(gamma);
    phi = R::rbeta(a_phi + sum_gamma, b_phi + p - sum_gamma);
    
    if((fmod(iter, thin) == 0.0) && (iter >= burn)){
      GAMMAVEC.row(iter - burn) = gamma.t();
      phivec(iter - burn) = phi;
      indvec(iter-burn) = ind;
    }
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("gamma") = GAMMAVEC,
    Rcpp::Named("prop") = indvec,
    Rcpp::Named("phi") = phivec,
    Rcpp::Named("loglik") = loglikvec,
    Rcpp::Named("acc_rate_bd") = acc_rate_bd,
    Rcpp::Named("acc_rate_move") = acc_rate_move,
    Rcpp::Named("acc_rate_bd_chain") = accratevec,
    Rcpp::Named("lambda_chain") = lambdavec
  );
  
}