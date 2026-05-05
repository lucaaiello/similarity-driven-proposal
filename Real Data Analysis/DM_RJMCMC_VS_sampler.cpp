// [[Rcpp::plugins(cpp14,openmp)]]
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <RcppEnsmallen.h>
#include <Rcpp.h>
#include <omp.h>
#include <sstream>
using namespace Rcpp;

// [[Rcpp::depends(RcppEnsmallen)]]
// [[Rcpp::depends(RcppArmadillo)]]

////////////////////////
// Optimization stuff //
////////////////////////

class DMRegressionFunction {
public:
  DMRegressionFunction(const arma::mat& X_j, 
                       const arma::imat& Y,
                       const arma::rowvec& intercepts,
                       const arma::ivec& PP,
                       const arma::mat& Xbeta,
                       const int j_aux,
                       const double lambda)
    : X_j(X_j), Y(Y), intercepts(intercepts), PP(PP),
      Xbeta(Xbeta), j_aux(j_aux),
      lambda(lambda), n(Y.n_rows), J(Y.n_cols)
  {
    
    P_j_aux = PP(j_aux);
    
  }
  
  //////////////////////////////////////
  // Objective and Gradient
  //////////////////////////////////////
  
  double EvaluateWithGradient(const arma::mat& beta_j, arma::mat& grad_flat) {
    // Handle empty block
    if (P_j_aux == 0) {
      grad_flat.zeros(0,1);
      return 0.0;
    }
    
    // Check input dimensions
    if ((int)beta_j.n_rows != P_j_aux || beta_j.n_cols != 1)
      Rcpp::stop("beta_j must be a column vector of length P_j_aux");
    
    grad_flat.zeros(P_j_aux,1);
    double loglik = 0.0;
    
    int n_threads = omp_get_max_threads();
    std::vector<arma::mat> grad_locals(n_threads, arma::mat(P_j_aux,1,arma::fill::zeros));
    std::vector<double> loglik_locals(n_threads, 0.0);
    std::vector<arma::rowvec> eta_locals(n_threads, arma::rowvec(J, arma::fill::zeros));
    std::vector<arma::rowvec> gamma_locals(n_threads, arma::rowvec(J, arma::fill::zeros));
    
    // Parallel loop over observations
#pragma omp parallel for schedule(static)
    for (int i = 0; i < (int)n; ++i) {
      int tid = omp_get_thread_num();
      arma::mat& grad_private = grad_locals[tid];
      double& loglik_private = loglik_locals[tid];
      arma::rowvec& eta = eta_locals[tid];
      arma::rowvec& gamma = gamma_locals[tid];
      
      const arma::irowvec& Yi = Y.row(i);
      
      double gamma_sum = 0.0;
      
      // Compute linear predictors and gamma
      for (int j = 0; j < J; ++j) {
        if (PP(j) == 0) {
          eta(j) = intercepts(j);
        } else if (j == j_aux) {
          eta(j) = intercepts(j_aux) + arma::dot(X_j.row(i), beta_j.col(0));
        } else {
          eta(j) = intercepts(j) + Xbeta(i, j); //arma::dot(X_list[j].row(i), beta_full[j]);
        }
        gamma(j) = std::exp(eta(j));
        gamma_sum += gamma(j);
      }
      
      int Ni = arma::accu(Yi);
      
      // Compute gradient only for j_aux
      double digamma_term = R::digamma(Yi(j_aux) + gamma(j_aux)) - R::digamma(gamma(j_aux))
        + R::digamma(gamma_sum) - R::digamma(gamma_sum + Ni);
      
      grad_private += -gamma(j_aux) * digamma_term * X_j.row(i).t();
      
      // Compute log-likelihood for observation i
      double ll_i = std::lgamma(Ni + 1) + std::lgamma(gamma_sum) - std::lgamma(Ni + gamma_sum);
      ll_i += arma::accu(arma::lgamma(Yi + gamma) - arma::lgamma(gamma) - arma::lgamma(Yi + 1));
      loglik_private += ll_i;
    }
    
    // Combine thread-local results
    for (int t = 0; t < n_threads; ++t) {
      grad_flat += grad_locals[t];
      loglik += loglik_locals[t];
    }
    
    // Apply ridge penalty
    grad_flat += lambda * beta_j;
    loglik -= 0.5 * lambda * arma::accu(beta_j % beta_j);
    
    // Return negative log-likelihood for minimization
    return -loglik;
  }
  
  
  //////////////////////////////////////
  // Hessian
  //////////////////////////////////////
  
  arma::mat ComputeHessian(const arma::mat& beta_j) {
    if (P_j_aux == 0)
      return arma::mat(0, 0, arma::fill::zeros);
    
    if (beta_j.n_cols != 1 || (int)beta_j.n_rows != P_j_aux)
      Rcpp::stop("beta_j must be a column vector of length P_j_aux");
    
    int n_threads = omp_get_max_threads();
    std::vector<arma::mat> H_locals(n_threads, arma::mat(P_j_aux, P_j_aux, arma::fill::zeros));
    std::vector<arma::rowvec> eta_locals(n_threads, arma::rowvec(J, arma::fill::zeros));
    std::vector<arma::rowvec> gamma_locals(n_threads, arma::rowvec(J, arma::fill::zeros));
    
    // Parallel loop over observations
#pragma omp parallel for schedule(static)
    for (int i = 0; i < (int)n; ++i) {
      int tid = omp_get_thread_num();
      arma::mat& H_local = H_locals[tid];
      arma::rowvec& eta = eta_locals[tid];
      arma::rowvec& gamma = gamma_locals[tid];
      
      const arma::irowvec& Yi = Y.row(i);
      
      // --- Compute η only for needed categories ---
      double Ai = 0.0;
      for (int j = 0; j < J; ++j) {
        if (PP(j) == 0) {
          eta(j) = intercepts(j);
          gamma(j) = std::exp(eta(j));
        } else if (j == j_aux) {
          eta(j) = intercepts(j_aux) + arma::dot(X_j.row(i), beta_j.col(0));
          gamma(j) = std::exp(eta(j));
        } else {
          // If beta not updated here, use current full beta
          eta(j) = intercepts(j) + Xbeta(i, j); // arma::dot(X_list[j].row(i), beta_full[j]);
          gamma(j) = std::exp(eta(j));
        }
        Ai += gamma(j);  // sum of all α for this observation
      }
      
      double Ni = arma::accu(Yi);
      
      // --- Gradient and second derivative for category j_aux ---
      double g_i = R::digamma(Yi(j_aux) + gamma(j_aux)) - R::digamma(gamma(j_aux))
        + R::digamma(Ai) - R::digamma(Ni + Ai);
      
      double t1 = R::trigamma(Yi(j_aux) + gamma(j_aux));
      double t2 = R::trigamma(gamma(j_aux));
      double tA = R::trigamma(Ai);
      double tNA = R::trigamma(Ni + Ai);
      
      double block_value = gamma(j_aux) * g_i + gamma(j_aux) * gamma(j_aux) * (t1 - t2 + tA - tNA);
      
      // --- Multiply by design matrix to get Hessian block ---
      const arma::rowvec& xj = X_j.row(i);  // no copy
      H_local += block_value * (xj.t() * xj);
    }
    
    // --- Combine thread-local contributions ---
    arma::mat H(P_j_aux, P_j_aux, arma::fill::zeros);
    for (int t = 0; t < n_threads; ++t)
      H += H_locals[t];
    
    // --- Apply ridge penalty and small numerical stabilization ---
    H -= lambda * arma::eye(P_j_aux, P_j_aux);
    
    H = arma::symmatu(H); // ensure symmetry
    
    return H;
  }
  
  
private:
  const arma::mat& X_j;
  const arma::imat& Y;
  const arma::rowvec& intercepts;
  const arma::ivec& PP;
  const arma::mat& Xbeta;
  const int j_aux;
  const double lambda;
  size_t n;
  int J;
  int P_j_aux;
};

//////////////////////////////////
// Wrappers for R
//////////////////////////////////

arma::vec dm_reg_lbfgs(const arma::mat& X_j, 
                       const arma::imat& Y,
                       const arma::rowvec& intercepts,
                       const arma::ivec& PP,
                       const arma::mat& Xbeta,
                       const int j,
                       const double lambda,
                       const arma::vec& beta_init,  
                       const int maxit = 10) { // 100
  
  DMRegressionFunction dmrf(X_j,Y,intercepts,PP,Xbeta,j,lambda);
  
  arma::vec beta_j = beta_init;
  
  ens::L_BFGS optimizer;
  optimizer.MaxIterations() = maxit;
  optimizer.Optimize(dmrf,beta_j);
  
  return beta_j;
}

arma::mat dm_reg_hessian(const arma::mat& X_j,
                         const arma::imat& Y,
                         const arma::rowvec& intercepts,
                         const arma::ivec& PP,
                         const arma::mat& Xbeta,
                         const int j,
                         const arma::vec& beta_j,
                         const double lambda) {
  
  DMRegressionFunction dmrf(X_j,Y,intercepts,PP,Xbeta,j,lambda);
  return dmrf.ComputeHessian(beta_j);
}


/////////////////////////
// Auxiliary variables //
/////////////////////////

// Convert binary vector to string key for hashing
std::string vectorToKey(const arma::ivec& config) {
  std::string key = "";
  for(arma::uword i = 0; i < config.n_elem; ++i) {
    key += std::to_string(config(i));
  }
  return key;
}

bool check_key_in_list(List mylist, std::string key) {
  CharacterVector names = mylist.names();
  
  for (int i = 0; i < names.size(); ++i) {
    if (as<std::string>(names[i]) == key) {
      return true;
    }
  }
  return false;
}

double dmvn_prec(const arma::vec& x,
                 const arma::vec& mean,
                 const arma::mat& Q,   // precision matrix
                 bool log_result = true) {
  int d = x.n_elem;
  arma::vec diff = x - mean;
  
  // Try Cholesky decomposition of Q (lower triangular)
  arma::mat L;
  bool chol_success = arma::chol(L, Q, "lower");
  
  double quadform, logdet;
  
  if (chol_success) {
    // Standard path: use Cholesky factor
    arma::vec y = L.t() * diff;
    quadform = arma::dot(y, y);
    logdet = 2.0 * arma::sum(arma::log(L.diag()));
  } else {
    // Fallback: add regularization and retry
    double reg = 1e-6;
    arma::mat Q_reg = Q;
    int max_attempts = 5;
    
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
      Q_reg.diag() += reg;
      chol_success = arma::chol(L, Q_reg, "lower");
      
      if (chol_success) {
        arma::vec y = L.t() * diff;
        quadform = arma::dot(y, y);
        logdet = 2.0 * arma::sum(arma::log(L.diag()));
        break;
      }
      reg *= 10.0;  // Increase regularization
    }
    
    if (!chol_success) {
      // Ultimate fallback: use diagonal approximation
      Rcpp::warning("Cholesky decomposition failed even with regularization. Using diagonal approximation.");
      arma::vec diag_Q = Q.diag();
      
      // Ensure positive diagonal elements
      diag_Q = arma::clamp(diag_Q, 1e-8, arma::datum::inf);
      
      quadform = arma::sum(diff % diff % diag_Q);
      logdet = arma::sum(arma::log(diag_Q));
    }
  }
  
  // log density
  double log_density = -0.5 * (d * std::log(2.0 * M_PI) - logdet + quadform);
  
  return log_result ? log_density : std::exp(log_density);
}

arma::vec rmvn_prec(const arma::vec& mu, const arma::mat& Q) {
  int d = mu.n_elem;
  
  // Try Cholesky decomposition of precision matrix Q = L * L^T
  arma::mat L;
  bool chol_success = arma::chol(L, Q, "lower");
  
  // Sample standard normal vector z ~ N(0, I)
  arma::vec z = arma::randn(d);
  
  if (chol_success) {
    // Standard path: solve L^T x = z
    arma::vec x = arma::solve(arma::trimatu(L.t()), z);
    return mu + x;
  } else {
    // Fallback: add regularization and retry
    double reg = 1e-6;
    arma::mat Q_reg = Q;
    int max_attempts = 5;
    
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
      Q_reg.diag() += reg;
      chol_success = arma::chol(L, Q_reg, "lower");
      
      if (chol_success) {
        arma::vec x = arma::solve(arma::trimatu(L.t()), z);
        return mu + x;
      }
      reg *= 10.0;  // Increase regularization
    }
    
    if (!chol_success) {
      // Ultimate fallback: use diagonal approximation
      Rcpp::warning("Cholesky decomposition failed even with regularization. Using diagonal approximation for sampling.");
      arma::vec diag_Q = Q.diag();
      
      // Ensure positive diagonal elements
      diag_Q = arma::clamp(diag_Q, 1e-8, arma::datum::inf);
      
      // Sample: x ~ N(mu, diag(1/Q_ii))
      arma::vec x = z / arma::sqrt(diag_Q);
      return mu + x;
    }
  }
  
  // Should never reach here
  return mu;
}

int sample_int(int n, 
               int size, 
               const arma::vec& probs, 
               bool replace = true) {
  int pn = probs.n_elem;
  
  if (pn != n) {
    Rcpp::stop("Length of probabilities vector must match n.");
  }
  // Convert arma::vec to Rcpp::NumericVector
  Rcpp::NumericVector rcpp_probs(probs.begin(), probs.end());
  
  // Perform sampling
  return Rcpp::sample(n, size, replace, rcpp_probs)[0];
}

double clamp_double(double value, double lower, double upper) {
  return std::max(lower, std::min(value, upper));
}

int sign_double(double value) {
  if (value > 0.0) {
    return 1;
  }
  if (value < 0.0) {
    return -1;
  }
  return 0;
}

arma::vec proposal_inf(const arma::ivec& xi_j, const int j,
                       const arma::ivec& PP, const arma::mat& Xbeta,
                       const std::unordered_map<std::string, arma::vec>& beta_j_hat_memory,
                       const arma::imat& Y, const arma::mat& X,
                       const arma::vec& alpha, const arma::mat& alpha_mat,
                       const arma::ivec& y_sum_i, 
                       const double log_lik_current,  
                       const arma::mat& gamma_current, 
                       double g, double c) {
  
  int P = xi_j.n_elem;
  int n = Y.n_rows;
  // double h = 0.01;
  
  arma::ivec Yj = Y.col(j);
  
  // ---- Baseline log-likelihood ----
  arma::mat gamma0 = arma::trunc_exp(alpha_mat);
  arma::vec gamma0_sum_i = arma::sum(gamma0, 1);
  
  double log_lik0 = arma::accu(
    arma::lgamma(gamma0_sum_i) - arma::lgamma(y_sum_i + gamma0_sum_i) +
      arma::sum(arma::lgamma(Y + gamma0) - arma::lgamma(gamma0), 1)
  );
  
  // Current state gamma sum
  arma::vec gamma_current_sum_i = arma::sum(gamma_current, 1);
  arma::vec gamma_j_current = gamma_current.col(j);
  arma::vec lgamma_gamma_current_sum_i = arma::lgamma(gamma_current_sum_i);
  arma::vec lgamma_y_sum_i_gamma_current_sum_i = arma::lgamma(y_sum_i + gamma_current_sum_i);
  arma::vec lgamma_Yj_gamma_j_current = arma::lgamma(Yj + gamma_j_current);
  arma::vec lgamma_gamma_j_current = arma::lgamma(gamma_j_current);
  
  // ---- Initialize proposal vectors ----
  arma::vec log_prop(P, arma::fill::zeros);
  
  // ---- Parallel loop ----
#pragma omp parallel for schedule(dynamic)
  for (int p = 0; p < P; p++) {
    
    arma::ivec xi_j_new = xi_j;
    xi_j_new(p) = 1 - xi_j(p);
    arma::uvec selected_new = arma::find(xi_j_new > 0);
    std::string config_new = vectorToKey(arma::conv_to<arma::ivec>::from(xi_j_new));
    
    arma::ivec PP_new = PP;
    PP_new(j) = static_cast<int>(selected_new.n_elem);
    double lambda_new = c * static_cast<double>(PP_new(j)) / static_cast<double>(n);
    
    // Subset the design matrix
    arma::mat X_new = X.cols(selected_new);
    
    // ---- Proposed state proposal (only if PP_prop(j) > 0) ----
    arma::vec beta_hat_j;
    
    if (PP_new(j) > 0) {
      auto it = beta_j_hat_memory.find(config_new);
      if (it != beta_j_hat_memory.end()) {
        // Warm start from previous beta_hat
        beta_hat_j = dm_reg_lbfgs(X_new, Y, alpha.t(), PP_new, Xbeta, j, lambda_new, it->second, 1);
      } else {
        // New configuration: start from random
        arma::vec beta_start = arma::randn(PP_new(j));
        beta_hat_j = dm_reg_lbfgs(X_new, Y, alpha.t(), PP_new, Xbeta, j, lambda_new, beta_start, 1);
      }
    } else {
      beta_hat_j = arma::vec();
    }
    
    // ---- Compute INCREMENTAL likelihood change (only for category j) ----
    arma::vec Xbeta_j_new(n);
    if (PP_new(j) > 0) {
      Xbeta_j_new = X_new * beta_hat_j;
    } else {
      Xbeta_j_new.zeros();
    }
    arma::vec gamma_j_new = arma::trunc_exp(alpha_mat.col(j) + Xbeta_j_new);
    
    // Change in row sums (current state -> proposed state)
    // arma::vec gamma_sum_i_new = gamma_current_sum_i - gamma_j_current + gamma_j_new;
    arma::vec gamma_sum_i_new = gamma_current_sum_i;
    gamma_sum_i_new += (gamma_j_new - gamma_j_current);
    
    // Incremental log-likelihood relative to CURRENT state
    // Δ1: change from category j terms
    arma::vec diff1 = arma::lgamma(Yj + gamma_j_new) - lgamma_Yj_gamma_j_current;
    arma::vec diff2 = lgamma_gamma_j_current - arma::lgamma(gamma_j_new);
    double Delta1 = arma::accu(diff1 + diff2);
    
    // Δ2: change from row sum terms
    arma::vec diff3 = arma::lgamma(gamma_sum_i_new) - lgamma_gamma_current_sum_i;
    arma::vec diff4 = lgamma_y_sum_i_gamma_current_sum_i - arma::lgamma(y_sum_i + gamma_sum_i_new);
    double Delta2 = arma::accu(diff3 + diff4);
    
    double log_lik_proposed = log_lik_current + Delta1 + Delta2;
    
    // ---- Likelihood ratio statistic (proposed vs intercept-only) ----
    double lik_ratio_stat = -2.0 * (log_lik0 - log_lik_proposed);
    if (!std::isfinite(lik_ratio_stat) || lik_ratio_stat < 0.0)
      lik_ratio_stat = 0.0;
    
    // ---- Compute -log10(p-value) ----
    double log_pval = R::pchisq(lik_ratio_stat, selected_new.n_elem, false, true); 
    double log10_LR_pvalue = -log_pval / std::log(10.0);
    
    // ---- Proposal weight ----
    log_prop(p) = arma::trunc_exp(g * arma::trunc_log(log10_LR_pvalue));
  }
  
  // ---- Normalize proposals ----
  arma::vec prop = arma::normalise(arma::trunc_exp(log_prop - log_prop.max()), 1);
  
  return prop;
}

double intercept_only_loglik(const arma::imat& Y,
                             const arma::mat& alpha_mat,
                             const arma::ivec& y_sum_i) {
  arma::mat gamma0 = arma::trunc_exp(alpha_mat);
  arma::vec gamma0_sum_i = arma::sum(gamma0, 1);
  
  return arma::accu(
    arma::lgamma(gamma0_sum_i) - arma::lgamma(y_sum_i + gamma0_sum_i) +
      arma::sum(arma::lgamma(Y + gamma0) - arma::lgamma(gamma0), 1)
  );
}

std::vector<int> active_with_inactive_neighbors(
    const arma::ivec& xi_j,
    const std::vector<std::vector<int>>& neigh_index) {
  std::vector<int> active_candidates;
  int P = xi_j.n_elem;
  
  for (int p = 0; p < P; ++p) {
    if (xi_j(p) != 1) {
      continue;
    }
    
    const std::vector<int>& neigh_p = neigh_index[p];
    bool has_inactive_neighbor = false;
    for (int q : neigh_p) {
      if (q >= 0 && q < P && xi_j(q) == 0) {
        has_inactive_neighbor = true;
        break;
      }
    }
    
    if (has_inactive_neighbor) {
      active_candidates.push_back(p);
    }
  }
  
  return active_candidates;
}

std::vector<int> inactive_neighbors_for_active(
    int active_idx,
    const arma::ivec& xi_j,
    const std::vector<std::vector<int>>& neigh_index) {
  std::vector<int> inactive_neighbors;
  int P = xi_j.n_elem;
  
  if (active_idx < 0 || active_idx >= P) {
    return inactive_neighbors;
  }
  
  const std::vector<int>& neigh_active = neigh_index[active_idx];
  inactive_neighbors.reserve(neigh_active.size());
  for (int q : neigh_active) {
    if (q >= 0 && q < P && xi_j(q) == 0) {
      inactive_neighbors.push_back(q);
    }
  }
  
  return inactive_neighbors;
}

arma::vec proposal_move_weights(
    const arma::ivec& xi_j,
    int drop_idx,
    const std::vector<int>& inactive_neighbors,
    int j,
    const arma::ivec& PP,
    const arma::mat& Xbeta,
    std::unordered_map<std::string, arma::vec>& beta_j_hat_memory,
    const arma::imat& Y,
    const arma::mat& X,
    const arma::vec& alpha,
    const arma::mat& alpha_mat,
    const arma::ivec& y_sum_i,
    double log_lik_current,
    const arma::mat& gamma_current,
    double lambda_move,
    double c) {
  int m = static_cast<int>(inactive_neighbors.size());
  if (m == 0) {
    return arma::vec();
  }
  
  int n = Y.n_rows;
  arma::ivec Yj = Y.col(j);
  double log_lik0 = intercept_only_loglik(Y, alpha_mat, y_sum_i);
  
  arma::vec gamma_current_sum_i = arma::sum(gamma_current, 1);
  arma::vec gamma_j_current = gamma_current.col(j);
  arma::vec lgamma_gamma_current_sum_i = arma::lgamma(gamma_current_sum_i);
  arma::vec lgamma_y_sum_i_gamma_current_sum_i =
    arma::lgamma(y_sum_i + gamma_current_sum_i);
  arma::vec lgamma_Yj_gamma_j_current = arma::lgamma(Yj + gamma_j_current);
  arma::vec lgamma_gamma_j_current = arma::lgamma(gamma_j_current);
  
  arma::vec score(m, arma::fill::zeros);
  
  for (int idx = 0; idx < m; ++idx) {
    int add_idx = inactive_neighbors[idx];
    arma::ivec xi_move = xi_j;
    xi_move(drop_idx) = 0;
    xi_move(add_idx) = 1;
    
    arma::uvec selected_move = arma::find(xi_move > 0);
    std::string config_move =
      vectorToKey(arma::conv_to<arma::ivec>::from(xi_move));
    
    arma::ivec PP_move = PP;
    PP_move(j) = static_cast<int>(selected_move.n_elem);
    double lambda_move_j =
      c * static_cast<double>(PP_move(j)) / static_cast<double>(n);
    
    arma::mat X_move = X.cols(selected_move);
    arma::vec beta_hat_move;
    
    if (PP_move(j) > 0) {
      auto it = beta_j_hat_memory.find(config_move);
      if (it != beta_j_hat_memory.end()) {
        beta_hat_move = dm_reg_lbfgs(
          X_move, Y, alpha.t(), PP_move, Xbeta, j, lambda_move_j, it->second, 1
        );
      } else {
        arma::vec beta_start = arma::randn(PP_move(j));
        beta_hat_move = dm_reg_lbfgs(
          X_move, Y, alpha.t(), PP_move, Xbeta, j, lambda_move_j, beta_start, 1
        );
      }
      beta_j_hat_memory[config_move] = beta_hat_move;
    } else {
      beta_hat_move = arma::vec();
    }
    
    arma::vec Xbeta_j_move(n, arma::fill::zeros);
    if (PP_move(j) > 0) {
      Xbeta_j_move = X_move * beta_hat_move;
    }
    
    arma::vec gamma_j_move = arma::trunc_exp(alpha_mat.col(j) + Xbeta_j_move);
    arma::vec gamma_sum_i_move = gamma_current_sum_i + gamma_j_move - gamma_j_current;
    
    arma::vec diff1 = arma::lgamma(Yj + gamma_j_move) - lgamma_Yj_gamma_j_current;
    arma::vec diff2 = lgamma_gamma_j_current - arma::lgamma(gamma_j_move);
    double delta1 = arma::accu(diff1 + diff2);
    
    arma::vec diff3 = arma::lgamma(gamma_sum_i_move) - lgamma_gamma_current_sum_i;
    arma::vec diff4 =
      lgamma_y_sum_i_gamma_current_sum_i - arma::lgamma(y_sum_i + gamma_sum_i_move);
    double delta2 = arma::accu(diff3 + diff4);
    
    double log_lik_move = log_lik_current + delta1 + delta2;
    double lik_ratio_stat = -2.0 * (log_lik0 - log_lik_move);
    if (!std::isfinite(lik_ratio_stat) || lik_ratio_stat < 0.0) {
      lik_ratio_stat = 0.0;
    }
    
    int df = std::max(1, static_cast<int>(selected_move.n_elem));
    double log_pval = R::pchisq(lik_ratio_stat, df, false, true);
    double log10_lrt_pvalue = -log_pval / std::log(10.0);
    score(idx) = arma::trunc_exp(lambda_move * arma::trunc_log(log10_lrt_pvalue));
  }
  
  return arma::normalise(arma::trunc_exp(score - score.max()), 1);
}

// [[Rcpp::export]]
Rcpp::List DM_informed_VS_RJMCMC(const arma::mat& X, const arma::imat& Y, 
                                 const arma::vec s2, const arma::vec r2,
                                 double a, double b, double g_initial, double c,
                                 int runs, int burn, int thin = 1, 
                                 bool informed = true,
                                 double trgt_acpt_alpha = 0.234,
                                 bool adapt_g = true,
                                 int adaptation_start = 100,
                                 int adaptation_end = -1,
                                 int window_size = 25,
                                 double rm_scale = 1.0,
                                 double decay_power = 0.75,
                                 double min_g = 0.05,
                                 double max_g = 10.0,
                                 Rcpp::Nullable<Rcpp::List> neigh = R_NilValue,
                                 bool use_move = false,
                                 int n_move_attempts = 1,
                                 double lambda_move = 1.0){
  
  if (runs <= 0) {
    Rcpp::stop("runs must be positive.");
  }
  if (burn < 0 || burn > runs) {
    Rcpp::stop("burn must lie between 0 and runs.");
  }
  if (thin < 0) {
    Rcpp::stop("thin must be non-negative.");
  }
  if (adapt_g) {
    if (window_size <= 0) {
      Rcpp::stop("window_size must be positive.");
    }
    if (rm_scale <= 0.0) {
      Rcpp::stop("rm_scale must be positive.");
    }
    if (decay_power <= 0.5 || decay_power > 1.0) {
      Rcpp::stop("decay_power must lie in (0.5, 1].");
    }
    if (min_g <= 0.0) {
      Rcpp::stop("min_g must be positive.");
    }
    if (max_g < min_g) {
      Rcpp::stop("max_g must be at least min_g.");
    }
  }
  if (use_move) {
    if (neigh.isNull()) {
      Rcpp::stop("neigh must be provided when use_move = true.");
    }
    if (n_move_attempts <= 0) {
      Rcpp::stop("n_move_attempts must be positive.");
    }
    if (lambda_move <= 0.0) {
      Rcpp::stop("lambda_move must be positive.");
    }
  }
  
  int adapt_start_iter = std::max(1, adaptation_start);
  int adapt_end_iter = (adaptation_end < 0) ? burn : adaptation_end;
  if (adapt_g) {
    adapt_end_iter = std::min(runs, adapt_end_iter);
    if (adapt_end_iter < adapt_start_iter) {
      Rcpp::stop("adaptation_end must be at least adaptation_start after truncation to the run length.");
    }
  }
  
  // const arma::vec& alpha,
  // Problem dimensions
  int n = Y.n_rows;
  int J = Y.n_cols;
  int P = X.n_cols;
  
  std::vector<std::vector<int>> neigh_index;
  if (use_move) {
    Rcpp::List neigh_list(neigh);
    if (neigh_list.size() != static_cast<R_xlen_t>(P)) {
      Rcpp::stop("neigh must have length equal to the number of predictors.");
    }
    
    neigh_index.resize(P);
    for (int p = 0; p < P; ++p) {
      arma::ivec neigh_p = Rcpp::as<arma::ivec>(neigh_list[p]);
      for (arma::uword idx = 0; idx < neigh_p.n_elem; ++idx) {
        int q = neigh_p(idx) - 1;
        if (q >= 0 && q < P && q != p) {
          neigh_index[p].push_back(q);
        }
      }
    }
    
    // Treat the neighborhood graph as undirected so forward and reverse
    // swap moves live on the same local structure.
    for (int p = 0; p < P; ++p) {
      std::vector<int> neigh_copy = neigh_index[p];
      for (int q : neigh_copy) {
        if (q >= 0 && q < P && q != p) {
          neigh_index[q].push_back(p);
        }
      }
    }
    
    for (int p = 0; p < P; ++p) {
      std::vector<int>& neigh_p = neigh_index[p];
      std::sort(neigh_p.begin(), neigh_p.end());
      neigh_p.erase(std::unique(neigh_p.begin(), neigh_p.end()), neigh_p.end());
    }
  }
  
  arma::ivec y_sum_i = arma::sum(Y, 1);
  
  // Unknown variables
  arma::vec alpha(J, arma::fill::zeros);
  arma::imat xi(P , J, arma::fill::zeros);
  arma::ivec PP = arma::sum(xi, 0).t();
  std::vector<std::string> config(J);
  
  std::vector<arma::uvec> selected(J);
  std::vector<arma::mat> X_list(J);
  std::vector<arma::vec> beta(J);
  arma::mat Xbeta(n, J);
  
  arma::vec log_prior_beta(J, arma::fill::zeros);
  arma::vec log_prior_xi(J, arma::fill::zeros);
  arma::vec lambda(J, arma::fill::zeros);
  
  std::vector<std::unordered_map<std::string, arma::vec>> beta_hat_memory(J);
  
  for (int j = 0; j < J; j++){
    config[j] = vectorToKey(arma::conv_to<arma::ivec>::from(xi.col(j)));
    selected[j] = arma::find(xi.col(j) > 0);
    X_list[j] = X.cols(selected[j]);
    beta[j] = arma::zeros(PP(j));    
    Xbeta.col(j) = X_list[j] * beta[j];
    log_prior_beta(j) = -(PP(j)/2.0) * (std::log(2.0 * M_PI) + arma::trunc_log(r2(j))) -
      0.5 * arma::dot(beta[j], beta[j])/r2(j);
    log_prior_xi(j) = arma::sum(arma::lgamma(xi.col(j) + a) + 
      arma::lgamma(1.0 - xi.col(j) + b)); 
    lambda(j) = c * static_cast<double>(PP(j)) / static_cast<double>(n);
  }
  
  arma::mat alpha_mat(n, J);
  alpha_mat.each_row() = alpha.t();
  arma::mat eta = alpha_mat + Xbeta;
  arma::mat gamma = arma::trunc_exp(eta);
  arma::vec gamma_sum_i = arma::sum(gamma, 1);
  
  // Compute matrix-wise lgamma terms
  arma::mat lgamma_Y_plus_gamma = arma::lgamma(Y + gamma);
  arma::mat lgamma_gamma = arma::lgamma(gamma);
  arma::vec term1 = arma::sum(lgamma_Y_plus_gamma - lgamma_gamma, 1);
  arma::vec term2 = arma::lgamma(gamma_sum_i) - arma::lgamma(y_sum_i + gamma_sum_i);
  
  // Total log-likelihood
  double log_lik = arma::accu(term1 + term2);
  
  double log_prior_alpha = -0.5 * arma::accu(arma::square(alpha) / s2);
  
  arma::mat prop_unif(P, J, arma::fill::ones);
  for (int j = 0; j < J; j++){
    prop_unif.col(j) /= P;
  }
  
  int post_runs;
  if(thin == 0){
    post_runs = floor(runs-burn);
  }else{
    post_runs = floor((runs-burn)/thin);
  }
  
  arma::mat alpha_post(post_runs, J, arma::fill::zeros);
  Rcpp::List beta_post(post_runs);
  arma::icube xi_post(P, J, post_runs, arma::fill::zeros);
  
  double s_alpha = 0.0;
  
  arma::vec m_alpha = alpha;
  
  arma::mat R_alpha = arma::eye(J, J); // 0.0005 *
  
  arma::mat Sigma_alpha = std::exp(s_alpha) * R_alpha;
  
  int accept_alpha = 0;
  arma::ivec accept_betaxi(J, arma::fill::zeros);
  int accept_betaxi_global = 0;
  arma::ivec accept_move(J, arma::fill::zeros);
  arma::ivec prop_move(J, arma::fill::zeros);
  int accept_move_global = 0;
  int prop_move_global = 0;
  
  // Category-wise Robbins-Monro state for g(j)
  arma::vec log_g(J, arma::fill::zeros);
  arma::vec prev_log_g(J, arma::fill::zeros);
  arma::vec prev_prev_log_g(J, arma::fill::zeros);
  arma::vec g(J, arma::fill::zeros);
  arma::vec prev_accept_rate_g(J, arma::fill::zeros);
  arma::ivec acc_betaxi_window(J, arma::fill::zeros);
  arma::ivec iter_counter_g(J, arma::fill::zeros);
  arma::ivec update_counter_g(J, arma::fill::zeros);
  int max_adapt_epochs = 0;
  if (adapt_g) {
    max_adapt_epochs = (adapt_end_iter - adapt_start_iter + 1) / window_size;
  }
  arma::mat g_window_acceptance(max_adapt_epochs, J);
  arma::mat g_window_values(max_adapt_epochs, J);
  g_window_acceptance.fill(arma::datum::nan);
  g_window_values.fill(arma::datum::nan);
  
  double log_g_min = std::log(min_g);
  double log_g_max = std::log(max_g);
  for (int j = 0; j < J; j++) {
    g(j) = clamp_double(g_initial, min_g, max_g);
    log_g(j) = std::log(g(j));
    prev_log_g(j) = log_g(j);
    prev_prev_log_g(j) = log_g(j);
  }
  
  // omp_set_num_threads(omp_get_max_threads() - 2);
  omp_set_num_threads(10);
  
#pragma omp parallel
{
#pragma omp single
{
  Rcpp::Rcout << "### Using " << omp_get_num_threads() << " threads \n" << std::endl;
}
}

auto chain_start = std::chrono::steady_clock::now();
int next_progress_hundredth_percent = 1;
std::size_t last_progress_width = 0;

for (int iter = 0; iter < runs; iter++) {
  int iter_id = iter + 1;
  
  double iter07 = pow(iter + 2, -0.7);
  
  //////////////////
  // alpha update //
  //////////////////
  
  arma::vec alpha_prop = arma::mvnrnd(alpha, Sigma_alpha);
  arma::mat alpha_mat_prop(n, J);
  alpha_mat_prop.each_row() = alpha_prop.t();
  
  arma::mat eta_prop = alpha_mat_prop + Xbeta;
  arma::mat gamma_prop = arma::trunc_exp(eta_prop);
  arma::vec gamma_prop_sum_i = arma::sum(gamma_prop, 1);
  
  double log_lik_prop = arma::accu(
    arma::lgamma(gamma_prop_sum_i) - arma::lgamma(y_sum_i + gamma_prop_sum_i) +
      arma::sum(arma::lgamma(Y + gamma_prop) - arma::lgamma(gamma_prop), 1)
  );
  
  double log_prior_alpha_prop = -0.5 * arma::accu(arma::square(alpha_prop) / s2);
  
  double log_ratio_alpha = log_lik_prop - log_lik + log_prior_alpha_prop - log_prior_alpha;
  
  double u = arma::randu();
  
  double log_acpt_alpha = std::min(0.0, log_ratio_alpha);
  
  if (log(u) < log_acpt_alpha) {
    alpha = alpha_prop;
    alpha_mat = alpha_mat_prop;
    eta = eta_prop;
    gamma = gamma_prop;
    gamma_sum_i = gamma_prop_sum_i;
    log_lik = log_lik_prop;
    log_prior_alpha = log_prior_alpha_prop;
    accept_alpha++;
  }
  
  // Update s_alpha
  s_alpha += iter07 * (std::exp(log_acpt_alpha) - trgt_acpt_alpha);
  
  // Calculate dalpha_2 and update m_alpha in a single loop
  arma::vec dalpha_2 = alpha - m_alpha;
  m_alpha += iter07 * dalpha_2;
  
  // Update R_alpha and xi_alpha using vectorized operations
  R_alpha = (1 - iter07) * R_alpha + iter07 * (dalpha_2 * dalpha_2.t());
  Sigma_alpha = std::exp(s_alpha) * R_alpha;
  
  // Adjust diagonal of xi_alpha using vectorized operations
  Sigma_alpha.diag() += 1e-3 * Sigma_alpha.diag();
  
  ///////////////////////
  // (beta, xi) update //
  ///////////////////////
  
  arma::vec log_prior_xi_prop(J, arma::fill::zeros);
  arma::vec log_prior_beta_prop(J, arma::fill::zeros);
  arma::vec lambda_prop(J, arma::fill::zeros);
  
  arma::imat xi_prop(P, J);
  arma::ivec PP_prop(J);
  std::vector<std::string> config_prop(J);
  std::vector<arma::uvec> selected_prop(J);
  std::vector<arma::mat> X_list_prop(J);
  std::vector<arma::vec> beta_prop(J);
  arma::mat Xbeta_prop(n, J);// = Xbeta;
  
  // Flag to track if any beta changed in this iteration
  bool any_betaxi_changed = false;
  
  for (int j = 0; j < J; j++){
    
    arma::vec proposal_xi_j(P);
    
    if (informed == true){
      proposal_xi_j = proposal_inf(xi.col(j), j, PP, Xbeta, beta_hat_memory[j], Y, X, alpha, alpha_mat, y_sum_i, log_lik, gamma, g(j), c);  
    } else {
      proposal_xi_j = prop_unif.col(j);  
    }
    
    int ind_j = sample_int(P, 1, proposal_xi_j) - 1;
    xi_prop.col(j) = xi.col(j);
    xi_prop(ind_j, j) = 1 - xi(ind_j, j);
    
    config_prop[j] = vectorToKey(arma::conv_to<arma::ivec>::from(xi_prop.col(j)));
    
    PP_prop = PP;
    PP_prop(j) = arma::sum(xi_prop.col(j));
    
    selected_prop[j] = arma::find(xi_prop.col(j) > 0);
    
    X_list_prop[j] = X.cols(selected_prop[j]);
    
    log_prior_xi_prop(j) = arma::sum(arma::lgamma(xi_prop.col(j) + a) + 
      arma::lgamma(1.0 - xi_prop.col(j) + b)); 
    
    lambda_prop(j) = c * static_cast<double>(PP_prop(j)) / static_cast<double>(n);
    
    // Proposed state proposal (only if PP_prop(j) > 0)
    arma::vec beta_j_prop_hat;
    arma::mat H_beta_j_prop;
    double log_proposal_beta_forward;
    
    Xbeta_prop = Xbeta;
    
    if (PP_prop(j) > 0) {
      auto it = beta_hat_memory[j].find(config_prop[j]);
      if (it != beta_hat_memory[j].end()) {
        // Warm start from previous beta_hat
        beta_j_prop_hat = dm_reg_lbfgs(X_list_prop[j], Y, alpha.t(), PP_prop, Xbeta, j, lambda_prop(j), it->second, 3);
      } else {
        // New configuration: start from random
        arma::vec beta_start = arma::randn(PP_prop(j));
        beta_j_prop_hat = dm_reg_lbfgs(X_list_prop[j], Y, alpha.t(), PP_prop, Xbeta, j, lambda_prop(j), beta_start, 3);
      }
      // Store optimized beta_hat in memory
      beta_hat_memory[j][config_prop[j]] = beta_j_prop_hat;
      
      arma::mat H_beta_j_prop = dm_reg_hessian(X_list_prop[j], Y, alpha.t(), PP_prop, Xbeta, j, beta_j_prop_hat, lambda_prop(j));
      beta_prop[j] = rmvn_prec(beta_j_prop_hat, -H_beta_j_prop);
      Xbeta_prop.col(j) = X_list_prop[j] * beta_prop[j];
      
      log_proposal_beta_forward = dmvn_prec(beta_prop[j], beta_j_prop_hat, -H_beta_j_prop);
    } else {
      beta_prop[j] = arma::vec();
      Xbeta_prop.col(j).zeros();
      log_proposal_beta_forward = 0.0;
    }
    
    // Reverse proposal (only if PP(j) > 0)
    double log_proposal_beta_reverse;
    arma::vec beta_j_hat;
    
    if (PP(j) > 0) {
      // Use stored beta_hat for the current configuration
      arma::vec beta_j_hat_reverse = dm_reg_lbfgs(X_list[j], Y, alpha.t(), PP, Xbeta_prop, j, lambda(j), beta_hat_memory[j][config[j]], 3);
      arma::mat H_beta_j_reverse = dm_reg_hessian(X_list[j], Y, alpha.t(), PP, Xbeta_prop, j, beta_j_hat_reverse, lambda(j));
      log_proposal_beta_reverse = dmvn_prec(beta[j], beta_j_hat_reverse, -H_beta_j_reverse);
    } else {
      log_proposal_beta_reverse = 0.0;
    }
    
    log_prior_beta_prop(j) = (PP_prop(j) > 0) ? 
    -(PP_prop(j)/2.0) * (std::log(2.0 * M_PI) + arma::trunc_log(r2(j))) - 
    0.5 * arma::dot(beta_prop[j], beta_prop[j])/r2(j) : 0.0;
    
    eta_prop.col(j) = alpha_mat.col(j) + Xbeta_prop.col(j);
    gamma_prop.col(j) = arma::trunc_exp(eta_prop.col(j));
    // ----- incremental gamma_prop_sum_i update -----
    gamma_prop_sum_i = gamma_sum_i + gamma_prop.col(j) - gamma.col(j); 
    
    // ----- incremental log-likelihood update (Δ1 + Δ2) -----
    // Δ1: change from column j contribution: sum_i [ lgamma(y_ij + g'_ij) - lgamma(g'_ij) - (lgamma(y_ij + g_ij) - lgamma(g_ij)) ]
    arma::ivec Yj_d = Y.col(j); // n x 1 as double
    arma::vec t_old = arma::lgamma(Yj_d + gamma.col(j)) - arma::lgamma(gamma.col(j));
    arma::vec t_new = arma::lgamma(Yj_d + gamma_prop.col(j)) - arma::lgamma(gamma_prop.col(j));
    double Delta1 = arma::accu(t_new - t_old);
    
    // Δ2: change from total-sum terms: sum_i [ lgamma(gsum') - lgamma(gsum) - lgamma(y_sum + gsum') + lgamma(y_sum + gsum) ]
    arma::vec s_old = arma::lgamma(gamma_sum_i) - arma::lgamma(y_sum_i + gamma_sum_i);
    arma::vec s_new = arma::lgamma(gamma_prop_sum_i) - arma::lgamma(y_sum_i + gamma_prop_sum_i);
    double Delta2 = arma::accu(s_new - s_old);
    
    double log_lik_prop = log_lik + Delta1 + Delta2;
    
    u = arma::randu();
    
    arma::vec proposal_xi_j_prop(P);
    
    if (informed == true) {
      proposal_xi_j_prop = proposal_inf(xi_prop.col(j), j, PP_prop, Xbeta_prop, beta_hat_memory[j], Y, X, alpha, alpha_mat, y_sum_i, log_lik_prop, gamma_prop, g(j), c);
    } else {
      proposal_xi_j_prop = prop_unif.col(j);  
    }
    
    double log_ratio_betaxi_j = log_lik_prop - log_lik +
      log_prior_beta_prop(j) - log_prior_beta(j) +
      log_prior_xi_prop(j) - log_prior_xi(j) +
      log_proposal_beta_reverse - log_proposal_beta_forward +
      arma::trunc_log(proposal_xi_j_prop(ind_j)) - arma::trunc_log(proposal_xi_j(ind_j));
    
    // double log_acpt_betaxi_j = std::min(0.0, log_ratio_betaxi_j);
    bool accepted_betaxi_j = false;
    
    if (log(u) < log_ratio_betaxi_j) {
      beta[j] = beta_prop[j];
      xi.col(j) = xi_prop.col(j);
      selected[j] = selected_prop[j];
      config[j] = config_prop[j];
      PP(j) = PP_prop(j);
      lambda(j) = lambda_prop(j);
      X_list[j] = X_list_prop[j];
      Xbeta.col(j) = Xbeta_prop.col(j);
      eta.col(j) = eta_prop.col(j);
      gamma.col(j) = gamma_prop.col(j);
      gamma_sum_i = gamma_prop_sum_i;
      log_lik = log_lik_prop;
      log_prior_beta(j) = log_prior_beta_prop(j);
      log_prior_xi(j) = log_prior_xi_prop(j);
      accept_betaxi(j)++;
      accepted_betaxi_j = true;
      // Mark that at least one beta changed
      any_betaxi_changed = true;
    }
    
    if (use_move && PP(j) > 0) {
      for (int move_attempt = 0; move_attempt < n_move_attempts; ++move_attempt) {
        std::vector<int> active_candidates =
          active_with_inactive_neighbors(xi.col(j), neigh_index);
        if (active_candidates.empty()) {
          break;
        }
        
        int active_pick = active_candidates[
        arma::randi(arma::distr_param(0, static_cast<int>(active_candidates.size()) - 1))
        ];
        std::vector<int> inactive_candidates =
          inactive_neighbors_for_active(active_pick, xi.col(j), neigh_index);
        if (inactive_candidates.empty()) {
          continue;
        }
        
        arma::vec move_weights = proposal_move_weights(
          xi.col(j), active_pick, inactive_candidates, j, PP, Xbeta,
          beta_hat_memory[j], Y, X, alpha, alpha_mat, y_sum_i,
          log_lik, gamma, lambda_move, c
        );
        if (move_weights.n_elem == 0) {
          continue;
        }
        
        int add_pos = sample_int(
          static_cast<int>(inactive_candidates.size()), 1, move_weights
        ) - 1;
        int add_idx = inactive_candidates[add_pos];
        
        prop_move(j) += 1;
        prop_move_global += 1;
        
        arma::ivec xi_move_j = xi.col(j);
        xi_move_j(active_pick) = 0;
        xi_move_j(add_idx) = 1;
        
        std::string config_move =
          vectorToKey(arma::conv_to<arma::ivec>::from(xi_move_j));
        arma::ivec PP_move = PP;
        PP_move(j) = arma::sum(xi_move_j);
        arma::uvec selected_move = arma::find(xi_move_j > 0);
        arma::mat X_move = X.cols(selected_move);
        
        double log_prior_xi_move = arma::sum(
          arma::lgamma(xi_move_j + a) + arma::lgamma(1.0 - xi_move_j + b)
        );
        double lambda_move_j =
          c * static_cast<double>(PP_move(j)) / static_cast<double>(n);
        
        arma::vec beta_j_move_hat;
        double log_proposal_beta_forward = 0.0;
        arma::vec beta_move_j;
        arma::mat Xbeta_move = Xbeta;
        
        if (PP_move(j) > 0) {
          auto it_move = beta_hat_memory[j].find(config_move);
          if (it_move != beta_hat_memory[j].end()) {
            beta_j_move_hat = dm_reg_lbfgs(
              X_move, Y, alpha.t(), PP_move, Xbeta, j, lambda_move_j,
              it_move->second, 3
            );
          } else {
            arma::vec beta_start = arma::randn(PP_move(j));
            beta_j_move_hat = dm_reg_lbfgs(
              X_move, Y, alpha.t(), PP_move, Xbeta, j, lambda_move_j,
              beta_start, 3
            );
          }
          beta_hat_memory[j][config_move] = beta_j_move_hat;
          
          arma::mat H_beta_j_move = dm_reg_hessian(
            X_move, Y, alpha.t(), PP_move, Xbeta, j, beta_j_move_hat, lambda_move_j
          );
          beta_move_j = rmvn_prec(beta_j_move_hat, -H_beta_j_move);
          Xbeta_move.col(j) = X_move * beta_move_j;
          log_proposal_beta_forward = dmvn_prec(
            beta_move_j, beta_j_move_hat, -H_beta_j_move
          );
        } else {
          beta_move_j = arma::vec();
          Xbeta_move.col(j).zeros();
        }
        
        double log_prior_beta_move = (PP_move(j) > 0) ?
        -(PP_move(j) / 2.0) * (std::log(2.0 * M_PI) + arma::trunc_log(r2(j))) -
        0.5 * arma::dot(beta_move_j, beta_move_j) / r2(j) :
          0.0;
        
        arma::vec eta_move_j = alpha_mat.col(j) + Xbeta_move.col(j);
        arma::vec gamma_move_j = arma::trunc_exp(eta_move_j);
        arma::vec gamma_sum_i_move = gamma_sum_i + gamma_move_j - gamma.col(j);
        
        arma::ivec Yj_move = Y.col(j);
        arma::vec t_old_move =
          arma::lgamma(Yj_move + gamma.col(j)) - arma::lgamma(gamma.col(j));
        arma::vec t_new_move =
          arma::lgamma(Yj_move + gamma_move_j) - arma::lgamma(gamma_move_j);
        double delta1_move = arma::accu(t_new_move - t_old_move);
        
        arma::vec s_old_move =
          arma::lgamma(gamma_sum_i) - arma::lgamma(y_sum_i + gamma_sum_i);
        arma::vec s_new_move =
          arma::lgamma(gamma_sum_i_move) - arma::lgamma(y_sum_i + gamma_sum_i_move);
        double delta2_move = arma::accu(s_new_move - s_old_move);
        
        double log_lik_move = log_lik + delta1_move + delta2_move;
        
        double log_proposal_beta_reverse = 0.0;
        if (PP(j) > 0) {
          auto it_reverse = beta_hat_memory[j].find(config[j]);
          arma::vec beta_start_reverse = (it_reverse != beta_hat_memory[j].end()) ?
          it_reverse->second : beta[j];
          arma::vec beta_j_hat_reverse = dm_reg_lbfgs(
            X_list[j], Y, alpha.t(), PP, Xbeta_move, j, lambda(j),
            beta_start_reverse, 3
          );
          arma::mat H_beta_j_reverse = dm_reg_hessian(
            X_list[j], Y, alpha.t(), PP, Xbeta_move, j, beta_j_hat_reverse, lambda(j)
          );
          log_proposal_beta_reverse = dmvn_prec(
            beta[j], beta_j_hat_reverse, -H_beta_j_reverse
          );
        }
        
        arma::mat gamma_move = gamma;
        gamma_move.col(j) = gamma_move_j;
        
        std::vector<int> active_reverse =
          active_with_inactive_neighbors(xi_move_j, neigh_index);
        auto active_reverse_it = std::find(
          active_reverse.begin(), active_reverse.end(), add_idx
        );
        if (active_reverse_it == active_reverse.end()) {
          continue;
        }
        
        std::vector<int> inactive_reverse =
          inactive_neighbors_for_active(add_idx, xi_move_j, neigh_index);
        auto inactive_reverse_it = std::find(
          inactive_reverse.begin(), inactive_reverse.end(), active_pick
        );
        if (inactive_reverse_it == inactive_reverse.end()) {
          continue;
        }
        
        arma::vec reverse_move_weights = proposal_move_weights(
          xi_move_j, add_idx, inactive_reverse, j, PP_move, Xbeta_move,
          beta_hat_memory[j], Y, X, alpha, alpha_mat, y_sum_i,
          log_lik_move, gamma_move, lambda_move, c
        );
        if (reverse_move_weights.n_elem == 0) {
          continue;
        }
        
        double pf = (1.0 / static_cast<double>(active_candidates.size())) *
          move_weights(add_pos);
        int reverse_add_pos = static_cast<int>(
          std::distance(inactive_reverse.begin(), inactive_reverse_it)
        );
        double pr = (1.0 / static_cast<double>(active_reverse.size())) *
          reverse_move_weights(reverse_add_pos);
        
        if (!std::isfinite(pf) || !std::isfinite(pr) || pf <= 0.0 || pr <= 0.0) {
          continue;
        }
        
        double log_ratio_move = log_lik_move - log_lik +
          log_prior_beta_move - log_prior_beta(j) +
          log_prior_xi_move - log_prior_xi(j) +
          log_proposal_beta_reverse - log_proposal_beta_forward +
          std::log(pr) - std::log(pf);
        
        if (std::log(arma::randu()) < log_ratio_move) {
          beta[j] = beta_move_j;
          xi.col(j) = xi_move_j;
          selected[j] = selected_move;
          config[j] = config_move;
          PP(j) = PP_move(j);
          lambda(j) = lambda_move_j;
          X_list[j] = X_move;
          Xbeta.col(j) = Xbeta_move.col(j);
          eta.col(j) = eta_move_j;
          gamma.col(j) = gamma_move_j;
          gamma_sum_i = gamma_sum_i_move;
          log_lik = log_lik_move;
          log_prior_beta(j) = log_prior_beta_move;
          log_prior_xi(j) = log_prior_xi_move;
          accept_move(j) += 1;
          accept_move_global += 1;
          break;
        }
      }
    }
    
    bool in_adaptation_window = adapt_g &&
      (iter_id >= adapt_start_iter) &&
      (iter_id <= adapt_end_iter);
    if (in_adaptation_window) {
      acc_betaxi_window(j) += (accepted_betaxi_j ? 1 : 0);
      iter_counter_g(j)++;
    }
    
  }
  
  // After the loop: increment counter if any beta changed
  if (any_betaxi_changed) {
    accept_betaxi_global++;  // Increment a single counter (not indexed by j)
  }
  
  ////////////////////////
  // ROBBINS-MONRO ADAPTATION FOR g (category-wise)
  ////////////////////////
  
  // This should be placed after the for loop over j (after all beta, xi updates)
  for (int j = 0; j < J; j++) {
    
    if (adapt_g && iter_counter_g(j) == window_size) {
      update_counter_g(j) += 1;
      double curr_accept_rate_g = static_cast<double>(acc_betaxi_window(j)) /
        static_cast<double>(window_size);
      double step_size = rm_scale *
        std::pow(static_cast<double>(update_counter_g(j)), -decay_power);
      
      // Bootstrap one bounded move after the first window, then use
      // the Robbins-Monro hill-climbing update thereafter.
      if (update_counter_g(j) == 1) {
        double probe_direction = (log_g(j) >= log_g_max) ? -1.0 : 1.0;
        double log_g_new = clamp_double(log_g(j) + probe_direction * step_size,
                                        log_g_min, log_g_max);
        prev_prev_log_g(j) = log_g(j);
        prev_log_g(j) = log_g_new;
        log_g(j) = log_g_new;
        
        // Acceptance change: Acc^(j-1) - Acc^(j-2)
        
        // Direction of last g change: I[log_g^(j-1) > log_g^(j-2)] - I[log_g^(j-1) ≤ log_g^(j-2)]
        
        // Gradient signal: positive if both increased together or both decreased together
        
        // Robbins-Monro update in log-space
        
        // Store previous for next iteration
        
        // Update with bounds
        
        // Convert back to regular space
        g(j) = std::exp(log_g(j));
        
      } else {
        double acc_change = curr_accept_rate_g - prev_accept_rate_g(j);
        double g_direction = static_cast<double>(
          sign_double(prev_log_g(j) - prev_prev_log_g(j))
        );
        double log_g_new = prev_log_g(j) + step_size * acc_change * g_direction;
        log_g_new = clamp_double(log_g_new, log_g_min, log_g_max);
        prev_prev_log_g(j) = prev_log_g(j);
        prev_log_g(j) = log_g_new;
        log_g(j) = log_g_new;
        g(j) = std::exp(log_g(j));
      }
      
      int epoch_idx = update_counter_g(j) - 1;
      if (epoch_idx >= 0 && epoch_idx < max_adapt_epochs) {
        g_window_acceptance(epoch_idx, j) = curr_accept_rate_g;
        g_window_values(epoch_idx, j) = g(j);
      }
      prev_accept_rate_g(j) = curr_accept_rate_g;
      acc_betaxi_window(j) = 0;
      iter_counter_g(j) = 0;
    }
  }
  
  // In the loop
  if((fmod(iter, thin) == 0.0) & (iter >= burn)){
    int save_idx = (iter - burn) / thin;  // Correct index for thinning
    alpha_post.row(save_idx) = alpha.t();
    beta_post[save_idx] = Rcpp::wrap(beta);  // Wrap converts std::vector<arma::vec> to R list
    xi_post.slice(save_idx) = xi;
  }
  
  double progress = 100.0 * static_cast<double>(iter + 1) /
    static_cast<double>(runs);
  int completed_progress_hundredth_percent =
    static_cast<int>(std::floor(progress * 100.0));
  double acc_rate_alpha = static_cast<double>(accept_alpha) /
    static_cast<double>(iter + 1);
  double acc_rate_betaxi_global = static_cast<double>(accept_betaxi_global) /
    static_cast<double>(iter + 1);
  double acc_rate_move_global = (prop_move_global > 0) ?
  static_cast<double>(accept_move_global) / static_cast<double>(prop_move_global) :
    0.0;
  arma::vec acc_rate_betaxi_by_category =
    arma::conv_to<arma::vec>::from(accept_betaxi) / static_cast<double>(iter + 1);
  double acc_rate_betaxi_cat_min = acc_rate_betaxi_by_category.min();
  double acc_rate_betaxi_cat_mean = arma::mean(acc_rate_betaxi_by_category);
  double acc_rate_betaxi_cat_max = acc_rate_betaxi_by_category.max();
  
  if (completed_progress_hundredth_percent >= next_progress_hundredth_percent) {
    double elapsed_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - chain_start).count();
    double average_iter_seconds = elapsed_seconds / static_cast<double>(iter + 1);
    double remaining_seconds = average_iter_seconds *
      static_cast<double>(runs - iter - 1);
    
    std::ostringstream progress_stream;
    progress_stream << std::fixed << std::setprecision(2)
                    << "Progress: " << progress
                    << "% | elapsed: " << elapsed_seconds
                    << "s | eta: " << remaining_seconds
                    << "s | acc_alpha: " << std::setprecision(3) << acc_rate_alpha
                    << " | acc_betaxi: " << acc_rate_betaxi_global
                    << " | acc_move: " << acc_rate_move_global
                    << " | cat[min/mean/max]: "
                    << acc_rate_betaxi_cat_min << "/"
                    << acc_rate_betaxi_cat_mean << "/"
                    << acc_rate_betaxi_cat_max;
    std::string progress_line = progress_stream.str();
    if (progress_line.size() < last_progress_width) {
      progress_line.append(last_progress_width - progress_line.size(), ' ');
    }
    last_progress_width = progress_line.size();
    
    Rprintf("\r%s", progress_line.c_str());
    R_FlushConsole();
    next_progress_hundredth_percent =
      completed_progress_hundredth_percent + 1;
  }
  
}

Rprintf("\n");
R_FlushConsole();

arma::vec acpt_rate_betaxi_by_category =
  arma::conv_to<arma::vec>::from(accept_betaxi) / static_cast<double>(runs);
arma::vec acpt_rate_move_by_category(J, arma::fill::zeros);
for (int j = 0; j < J; ++j) {
  if (prop_move(j) > 0) {
    acpt_rate_move_by_category(j) =
      static_cast<double>(accept_move(j)) / static_cast<double>(prop_move(j));
  }
}

return Rcpp::List::create(
  Rcpp::Named("alpha") = alpha_post,
  Rcpp::Named("beta") = beta_post,
  Rcpp::Named("xi") = xi_post,
  Rcpp::Named("acpt_rate_alpha") = static_cast<double>(accept_alpha) / static_cast<double>(runs),
  Rcpp::Named("acpt_rate_betaxi_global") = static_cast<double>(accept_betaxi_global) / static_cast<double>(runs),
  Rcpp::Named("acpt_rate_betaxi_by_category") = acpt_rate_betaxi_by_category,
  Rcpp::Named("acpt_count_betaxi_by_category") = accept_betaxi,
  Rcpp::Named("acpt_rate_move_global") = (prop_move_global > 0) ?
static_cast<double>(accept_move_global) / static_cast<double>(prop_move_global) : 0.0,
Rcpp::Named("acpt_rate_move_by_category") = acpt_rate_move_by_category,
Rcpp::Named("acpt_count_move_by_category") = accept_move,
Rcpp::Named("prop_count_move_by_category") = prop_move,
Rcpp::Named("g_final") = g,
Rcpp::Named("g_window_acceptance") = g_window_acceptance,
Rcpp::Named("g_window_values") = g_window_values
);

}