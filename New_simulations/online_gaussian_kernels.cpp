#include <RcppArmadillo.h>
#include "gaussian_revision_shared.h"
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

struct Counters {
  double target_seconds = 0.0;
  double auxiliary_seconds = 0.0;
  double forward_seconds = 0.0;
  double reverse_seconds = 0.0;
  long long target_evaluations = 0;
  long long ols_fits = 0;
  long long proposal_normalizations = 0;
};

static double elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

static arma::uvec selected_columns(const arma::uvec& gamma) {
  return arma::join_cols(arma::uvec({0}), arma::find(gamma > 0) + 1);
}

static double log_beta_binomial_configuration(const arma::uvec& gamma,
                                               double a_phi, double b_phi) {
  double k = arma::accu(gamma);
  double p = gamma.n_elem;
  return R::lbeta(a_phi + k, b_phi + p - k) - R::lbeta(a_phi, b_phi);
}

static double exact_log_target_raw(const arma::uvec& gamma,
                                   const arma::mat& X, const arma::vec& y,
                                   double a, double d, double tau,
                                   double a_phi, double b_phi) {
  arma::uvec selected = selected_columns(gamma);
  arma::mat Xm = X.cols(selected);
  arma::mat precision = Xm.t() * Xm + arma::eye(selected.n_elem, selected.n_elem) / tau;
  arma::vec rhs = Xm.t() * y;
  arma::mat L = arma::chol(precision, "lower");
  arma::vec mean = arma::solve(arma::trimatu(L.t()),
                    arma::solve(arma::trimatl(L), rhs));
  double astar = a + 0.5 * (arma::dot(y, y) - arma::dot(mean, rhs));
  double dstar = d + 0.5 * y.n_elem;
  double logdet = 2.0 * arma::sum(arma::log(L.diag()));
  double loglik = -0.5 * y.n_elem * (std::log(M_PI) + std::log(2.0))
    + 0.5 * (-selected.n_elem * std::log(tau) - logdet)
    + d * std::log(a) - dstar * std::log(astar)
    + std::lgamma(dstar) - std::lgamma(d);
  return loglik + log_beta_binomial_configuration(gamma, a_phi, b_phi);
}

static double exact_log_target(const arma::uvec& gamma,
                               const arma::mat& X, const arma::vec& y,
                               double a, double d, double tau,
                               double a_phi, double b_phi, Counters& counters) {
  auto start = std::chrono::steady_clock::now();
  double value = exact_log_target_raw(gamma, X, y, a, d, tau, a_phi, b_phi);
  counters.target_seconds += elapsed_seconds(start);
  counters.target_evaluations += 1;
  return value;
}

static double model_rss(const arma::uvec& gamma, const arma::mat& X,
                        const arma::vec& y, Counters& counters) {
  auto start = std::chrono::steady_clock::now();
  arma::uvec selected = selected_columns(gamma);
  arma::mat Xm = X.cols(selected);
  double rss = gaussrev::rss(Xm.t() * Xm, Xm.t() * y, arma::dot(y, y));
  counters.auxiliary_seconds += elapsed_seconds(start);
  counters.ols_fits += 1;
  return std::max(rss, 1e-14);
}

static double f_similarity_score(const arma::uvec& gamma, const arma::mat& X,
                                 const arma::vec& y, double rss_null,
                                 Counters& counters) {
  int k = arma::accu(gamma);
  if (k == 0) return 0.0;
  double rss = model_rss(gamma, X, y, counters);
  int df2 = static_cast<int>(y.n_elem) - k - 1;
  if (df2 <= 0) return 0.0;
  double f = std::max(((rss_null - rss) / k) / (rss / df2), 0.0);
  double logp = R::pf(f, 1, df2, false, true);
  return gaussrev::tail_score(logp);
}

static double bic_similarity_score(const arma::uvec& gamma, const arma::mat& X,
                                   const arma::vec& y, double rss_null,
                                   Counters& counters) {
  int k = arma::accu(gamma);
  if (k == 0) return 0.0;
  double rss = model_rss(gamma, X, y, counters);
  double lr = y.n_elem * (std::log(rss_null) - std::log(rss));
  double penalized = std::max(lr - std::log(static_cast<double>(y.n_elem)) * k, 0.0);
  double logp = R::pchisq(penalized, 1.0, false, true);
  return gaussrev::tail_score(logp);
}

static double bic_approximate_log_posterior(const arma::uvec& gamma,
                                            const arma::mat& X, const arma::vec& y,
                                            double a_phi, double b_phi,
                                            Counters& counters) {
  double rss = model_rss(gamma, X, y, counters);
  int parameters = arma::accu(gamma) + 1;
  double bic_part = -0.5 * (y.n_elem * std::log(rss / y.n_elem)
                            + parameters * std::log(static_cast<double>(y.n_elem)));
  return bic_part + log_beta_binomial_configuration(gamma, a_phi, b_phi);
}

static double log_barker(double delta) {
  if (delta > 0) return -std::log1p(std::exp(-delta));
  return delta - std::log1p(std::exp(delta));
}

static arma::vec proposal_probabilities(const arma::uvec& gamma,
                                        const std::string& method,
                                        const arma::mat& X, const arma::vec& y,
                                        double rss_null, double lambda,
                                        double current_log_target,
                                        double a, double d, double tau,
                                        double a_phi, double b_phi,
                                        Counters& counters) {
  int p = gamma.n_elem;
  arma::vec log_weights(p, arma::fill::zeros);
  double current_bic = 0.0;
  if (method == "bic_locally_balanced") {
    current_bic = bic_approximate_log_posterior(gamma, X, y, a_phi, b_phi, counters);
  }
  for (int j = 0; j < p; ++j) {
    arma::uvec candidate = gamma;
    candidate(j) = 1u - candidate(j);
    // Auxiliary-score limitations do not restrict posterior support.
    if (method == "uniform") {
      log_weights(j) = 0.0;
    } else if (method == "f_similarity") {
      log_weights(j) = std::pow(f_similarity_score(candidate, X, y, rss_null, counters), lambda);
    } else if (method == "bic_similarity") {
      log_weights(j) = std::pow(bic_similarity_score(candidate, X, y, rss_null, counters), lambda);
    } else if (method == "exact_locally_balanced") {
      double candidate_target = exact_log_target(candidate, X, y, a, d, tau, a_phi, b_phi, counters);
      log_weights(j) = log_barker(candidate_target - current_log_target);
    } else if (method == "bic_locally_balanced") {
      double candidate_bic = bic_approximate_log_posterior(candidate, X, y, a_phi, b_phi, counters);
      log_weights(j) = log_barker(candidate_bic - current_bic);
    } else {
      stop("Unknown method: " + method);
    }
  }
  double maximum = log_weights.max();
  arma::vec probabilities = arma::exp(log_weights - maximum);
  probabilities /= arma::accu(probabilities);
  counters.proposal_normalizations += 1;
  return probabilities;
}

// Compute every one-flip F/BIC score from one factorization of the current
// model.  Addition RSS values use partial-regression identities; deletion RSS
// values use the diagonal of the inverse Gram matrix.
static arma::vec fast_similarity_probabilities(const arma::uvec& gamma,
                                                const std::string& score_type,
                                                const arma::mat& gram,
                                                const arma::vec& cross_y,
                                                double yty, double rss_null,
                                                int n, double lambda,
                                                bool balance_directions,
                                                double uniform_mixture,
                                                Counters& counters) {
  auto start = std::chrono::steady_clock::now();
  int p = gamma.n_elem;
  int k = arma::accu(gamma);
  arma::uvec selected = selected_columns(gamma);
  arma::mat gram_inv;
  arma::vec beta;
  double rss_current = 0.0;
  arma::ivec active_position(p, arma::fill::value(-1));
  bool rank_one_ok = false;
  if (score_type != "marginal") {
    arma::mat gram_s = gram.submat(selected, selected);
    rank_one_ok = selected.n_elem < (arma::uword)n && arma::inv_sympd(gram_inv, gram_s);
    if (!rank_one_ok) gram_inv = arma::pinv(gram_s);
    beta = gram_inv * cross_y(selected);
    rss_current = std::max(yty - arma::dot(cross_y(selected), beta), 1e-14);
    for (arma::uword pos = 1; pos < selected.n_elem; ++pos)
      active_position(selected(pos) - 1) = static_cast<int>(pos);
  }

  arma::vec log_weights(p, arma::fill::zeros);
  arma::vec marginal_evidence;
  double active_marginal_sum = 0.0;
  if (score_type == "marginal") {
    marginal_evidence.set_size(p);
    double intercept_gram = gram(0,0);
    for (int j = 0; j < p; ++j) {
      int column = j + 1;
      double sxx = gram(column,column) - gram(0,column)*gram(0,column)/intercept_gram;
      double sxy = cross_y(column) - gram(0,column)*cross_y(0)/intercept_gram;
      double explained = sxx > 1e-12 ? sxy*sxy/sxx : 0.0;
      double residual = std::max(rss_null-explained,1e-14);
      double f = std::max(explained/(residual/std::max(n-2,1)),0.0);
      marginal_evidence(j) = std::max(-R::pf(f,1,std::max(n-2,1),false,true)/std::log(10.0),0.0);
      if (gamma(j)) active_marginal_sum += marginal_evidence(j);
    }
  }
  for (int j = 0; j < p; ++j) {
    int candidate_k = k + (gamma(j) ? -1 : 1);
    double candidate_rss = 0.0;
    if (score_type != "marginal") {
      if (!rank_one_ok) {
        arma::uvec candidate = gamma; candidate(j) = 1u - candidate(j);
        arma::uvec cols = selected_columns(candidate);
        candidate_rss = gaussrev::rss(gram.submat(cols, cols), cross_y(cols), yty);
      } else if (gamma(j)) {
        int pos = active_position(j);
        candidate_rss = rss_current + beta(pos) * beta(pos) / gram_inv(pos, pos);
      } else {
        int column = j + 1;
        arma::vec cross = gram.submat(selected, arma::uvec({static_cast<arma::uword>(column)}));
        double denominator = gram(column, column) - arma::as_scalar(cross.t() * gram_inv * cross);
        double numerator = cross_y(column) - arma::dot(cross, beta);
        candidate_rss = denominator > 1e-12
          ? rss_current - numerator * numerator / denominator : rss_current;
      }
      candidate_rss = std::max(candidate_rss, 1e-14);
    }
    double score = 0.0;
    if (candidate_k > 0 && score_type == "f") {
      score = gaussrev::f_score(candidate_k, n, rss_null, candidate_rss);
    } else if (candidate_k > 0 && score_type == "bic") {
      score = gaussrev::lr_score(candidate_k, n, rss_null, candidate_rss, "bic");
    } else if (candidate_k > 0 && score_type == "marginal") {
      double candidate_sum = active_marginal_sum + (gamma(j) ? -marginal_evidence(j) : marginal_evidence(j));
      // A deliberately weak screening score: aggregate univariate evidence,
      // normalized only by sqrt(model size).  It ignores conditional effects.
      score = candidate_sum/std::sqrt(static_cast<double>(candidate_k));
    } else if (score_type != "f" && score_type != "bic" && score_type != "marginal") {
      stop("score_type must be f, bic, or marginal");
    }
    log_weights(j) = std::pow(score, lambda);
  }
  arma::vec probabilities(p, arma::fill::zeros);
  if (!balance_directions) {
    double maximum = log_weights.max();
    probabilities = arma::exp(log_weights - maximum);
    probabilities /= arma::accu(probabilities);
  } else {
    std::vector<arma::uword> birth_storage, death_storage;
    for (int j = 0; j < p; ++j) {
      if (!std::isfinite(log_weights(j))) continue;
      if (gamma(j)) death_storage.push_back(j); else birth_storage.push_back(j);
    }
    arma::uvec births(birth_storage), deaths(death_storage);
    bool has_births = births.n_elem > 0, has_deaths = deaths.n_elem > 0;
    double birth_mass = has_births && has_deaths ? 0.5 : (has_births ? 1.0 : 0.0);
    double death_mass = has_births && has_deaths ? 0.5 : (has_deaths ? 1.0 : 0.0);
    if (has_births) {
      arma::vec z = log_weights.elem(births);
      z = arma::exp(z - z.max());
      probabilities.elem(births) = birth_mass * z / arma::accu(z);
    }
    if (has_deaths) {
      arma::vec z = log_weights.elem(deaths);
      z = arma::exp(z - z.max());
      probabilities.elem(deaths) = death_mass * z / arma::accu(z);
    }
  }
  if (uniform_mixture > 0.0) {
    if (uniform_mixture >= 1.0) stop("uniform_mixture must be smaller than one");
    if (balance_directions) {
      std::vector<arma::uword> birth_storage, death_storage;
      for (int j = 0; j < p; ++j) {
        if (!std::isfinite(log_weights(j))) continue;
        if (gamma(j)) death_storage.push_back(j); else birth_storage.push_back(j);
      }
      arma::uvec births(birth_storage), deaths(death_storage);
      bool has_births = births.n_elem > 0, has_deaths = deaths.n_elem > 0;
      double birth_mass = has_births && has_deaths ? 0.5 : (has_births ? 1.0 : 0.0);
      double death_mass = has_births && has_deaths ? 0.5 : (has_deaths ? 1.0 : 0.0);
      probabilities *= (1.0 - uniform_mixture);
      if (has_births) probabilities.elem(births) += birth_mass * uniform_mixture / births.n_elem;
      if (has_deaths) probabilities.elem(deaths) += death_mass * uniform_mixture / deaths.n_elem;
    } else {
      arma::uvec valid = arma::find_finite(log_weights);
      probabilities *= (1.0 - uniform_mixture);
      probabilities.elem(valid) += uniform_mixture / valid.n_elem;
    }
  }
  counters.auxiliary_seconds += elapsed_seconds(start);
  if (score_type != "marginal")
    counters.ols_fits += 1; // one current-model factorization, not p candidate fits
  counters.proposal_normalizations += 1;
  return probabilities;
}

// Rank-one implementation of the comparator proposals.  Exact locally
// balanced weights use determinant/Schur-complement identities for the
// conjugate marginal target; BIC locally balanced weights use the analogous
// OLS RSS identities.  Thus each full one-flip neighborhood requires one
// current-model factorization rather than p separate fits.
static arma::vec fast_comparator_probabilities(const arma::uvec& gamma,
                                               const std::string& method,
                                               const arma::mat& gram,
                                               const arma::vec& cross_y,
                                               double yty, int n,
                                               double a, double d, double tau,
                                               double a_phi, double b_phi,
                                               Counters& counters) {
  const int p = gamma.n_elem;
  if (method == "uniform") {
    counters.proposal_normalizations += 1;
    return arma::vec(p, arma::fill::value(1.0 / p));
  }
  if (method != "exact_locally_balanced" && method != "bic_locally_balanced")
    stop("Fast comparator must be uniform, exact_locally_balanced, or bic_locally_balanced");

  auto start = std::chrono::steady_clock::now();
  arma::uvec selected = selected_columns(gamma);
  int current_k = arma::accu(gamma);
  arma::ivec active_position(p, arma::fill::value(-1));
  for (arma::uword pos = 1; pos < selected.n_elem; ++pos)
    active_position(selected(pos) - 1) = static_cast<int>(pos);
  arma::vec log_weights(p, arma::fill::zeros);

  if (method == "exact_locally_balanced") {
    arma::mat precision = gram.submat(selected, selected)
      + arma::eye(selected.n_elem, selected.n_elem) / tau;
    arma::mat precision_inverse;
    if (!arma::inv_sympd(precision_inverse, precision))
      stop("Current-model posterior precision is not positive definite");
    double sign = 0.0, logdet = 0.0;
    arma::log_det(logdet, sign, precision);
    if (sign <= 0.0) stop("Current-model posterior precision has nonpositive determinant");
    arma::vec rhs = cross_y(selected);
    arma::vec mean = precision_inverse * rhs;
    double quad = arma::dot(rhs, mean);
    double dstar = d + 0.5 * n;
    auto target_from_parts = [&](int model_k, double candidate_logdet, double candidate_quad) {
      int q = model_k + 1;
      double astar = a + 0.5 * (yty - candidate_quad);
      double loglik = -0.5 * n * (std::log(M_PI) + std::log(2.0))
        + 0.5 * (-q * std::log(tau) - candidate_logdet)
        + d * std::log(a) - dstar * std::log(astar)
        + std::lgamma(dstar) - std::lgamma(d);
      double logprior = R::lbeta(a_phi + model_k, b_phi + p - model_k)
        - R::lbeta(a_phi, b_phi);
      return loglik + logprior;
    };
    double current_target = target_from_parts(current_k, logdet, quad);
    for (int j = 0; j < p; ++j) {
      int candidate_k = current_k + (gamma(j) ? -1 : 1);
      double candidate_logdet, candidate_quad;
      if (gamma(j)) {
        int pos = active_position(j);
        double inverse_diagonal = precision_inverse(pos, pos);
        candidate_logdet = logdet + std::log(inverse_diagonal);
        candidate_quad = quad - mean(pos) * mean(pos) / inverse_diagonal;
      } else {
        arma::uword column = static_cast<arma::uword>(j + 1);
        arma::vec cross = gram.submat(selected, arma::uvec({column}));
        double schur = gram(column, column) + 1.0 / tau
          - arma::as_scalar(cross.t() * precision_inverse * cross);
        if (schur <= 1e-12) stop("Nonpositive Schur complement in exact locally balanced proposal");
        double residual = cross_y(column) - arma::dot(cross, mean);
        candidate_logdet = logdet + std::log(schur);
        candidate_quad = quad + residual * residual / schur;
      }
      double candidate_target = target_from_parts(candidate_k, candidate_logdet, candidate_quad);
      log_weights(j) = log_barker(candidate_target - current_target);
    }
    counters.target_evaluations += p + 1;
    counters.ols_fits += 1;
  } else {
    arma::mat gram_current = gram.submat(selected, selected);
    arma::mat gram_inverse;
    bool rank_one_ok = selected.n_elem < (arma::uword)n && arma::inv_sympd(gram_inverse, gram_current);
    if (!rank_one_ok) gram_inverse = arma::pinv(gram_current);
    arma::vec rhs = cross_y(selected);
    arma::vec beta = gram_inverse * rhs;
    double rss_current = std::max(yty - arma::dot(rhs, beta), 1e-14);
    auto bic_target = [&](int model_k, double rss) {
      int parameters = model_k + 1;
      return -0.5 * (n * std::log(rss / n) + parameters * std::log(static_cast<double>(n)))
        + R::lbeta(a_phi + model_k, b_phi + p - model_k) - R::lbeta(a_phi, b_phi);
    };
    double current_bic = bic_target(current_k, rss_current);
    for (int j = 0; j < p; ++j) {
      int candidate_k = current_k + (gamma(j) ? -1 : 1);
      double candidate_rss;
      if (!rank_one_ok) {
        arma::uvec candidate = gamma; candidate(j) = 1u - candidate(j);
        arma::uvec cols = selected_columns(candidate);
        candidate_rss = gaussrev::rss(gram.submat(cols, cols), cross_y(cols), yty);
      } else if (gamma(j)) {
        int pos = active_position(j);
        candidate_rss = rss_current + beta(pos) * beta(pos) / gram_inverse(pos, pos);
      } else {
        arma::uword column = static_cast<arma::uword>(j + 1);
        arma::vec cross = gram.submat(selected, arma::uvec({column}));
        double denominator = gram(column, column)
          - arma::as_scalar(cross.t() * gram_inverse * cross);
        double residual = cross_y(column) - arma::dot(cross, beta);
        candidate_rss = denominator > 1e-12
          ? rss_current - residual * residual / denominator : rss_current;
      }
      candidate_rss = std::max(candidate_rss, 1e-14);
      log_weights(j) = log_barker(bic_target(candidate_k, candidate_rss) - current_bic);
    }
    counters.ols_fits += 1;
  }
  double maximum = log_weights.max();
  arma::vec probabilities = arma::exp(log_weights - maximum);
  probabilities /= arma::accu(probabilities);
  counters.auxiliary_seconds += elapsed_seconds(start);
  counters.proposal_normalizations += 1;
  return probabilities;
}

// [[Rcpp::export]]
Rcpp::List run_optimized_similarity_scaling(const arma::vec& y,
                                             const arma::mat& predictors,
                                             std::string score_type,
                                             arma::uvec gamma,
                                             double time_budget_seconds,
                                             int burn = 500,
                                             int max_iterations = 1000000,
                                             double lambda_initial = 0.7,
                                             bool adapt_lambda = true,
                                             double target_acceptance = 0.234,
                                             int adaptation_window = 25,
                                             double lambda_min = 0.01,
                                             double lambda_max = 2.0,
                                             double a = 2.0, double d = 2.0,
                                             double tau = 1.0,
                                             double a_phi = 1.0, double b_phi = 3.0,
                                             bool balance_directions = false,
                                             std::string adaptation_metric = "acceptance",
                                             double target_effective_neighbors = -1.0,
                                             double adaptation_gain = 0.5,
                                             double adaptation_decay = 0.6,
                                             double adaptation_offset = 100.0,
                                             double uniform_mixture = 0.0) {
  RNGScope scope;
  if (score_type != "f" && score_type != "bic" && score_type != "marginal")
    stop("score_type must be f, bic, or marginal");
  if (time_budget_seconds <= 0 || burn < 0 || max_iterations <= burn)
    stop("Invalid time budget, burn, or maximum iteration count");
  if (lambda_initial <= 0 || lambda_min <= 0 || lambda_max < lambda_min)
    stop("Invalid lambda values");
  if (adaptation_window < 1) stop("adaptation_window must be positive");
  if (adaptation_metric != "acceptance" && adaptation_metric != "diversity")
    stop("adaptation_metric must be acceptance or diversity");
  if (adaptation_decay <= 0.5 || adaptation_decay > 1.0)
    stop("adaptation_decay must be in (0.5, 1]");
  if (uniform_mixture < 0.0 || uniform_mixture >= 1.0)
    stop("uniform_mixture must be in [0, 1)");
  if (target_effective_neighbors <= 0.0)
    target_effective_neighbors = std::max(5.0, 0.025 * gamma.n_elem);
  auto total_start = std::chrono::steady_clock::now();
  arma::mat X = arma::join_rows(arma::ones(y.n_elem), predictors);
  arma::mat gram = X.t() * X;
  arma::vec cross_y = X.t() * y;
  double yty = arma::dot(y, y);
  double rss_null = arma::accu(arma::square(y - arma::mean(y)));
  Counters counters;
  arma::vec pip_sum(gamma.n_elem, arma::fill::zeros);
  std::vector<double> model_size_trace, log_target_trace, lambda_trace, model_index_trace;
  int accepted = 0, window_accepted = 0, window_attempted = 0, updates = 0, iterations = 0;
  double window_proposal_ess_sum = 0.0;
  int retained_accepted = 0, birth_attempts = 0, birth_accepted = 0;
  int death_attempts = 0, death_accepted = 0;
  double proposal_ess_sum = 0.0, proposal_entropy_sum = 0.0, proposal_max_sum = 0.0;
  double log_lambda = std::log(lambda_initial);
  double current_target = exact_log_target(gamma, X, y, a, d, tau, a_phi, b_phi, counters);
  arma::vec current_proposal = fast_similarity_probabilities(
    gamma, score_type, gram, cross_y, yty, rss_null, y.n_elem, std::exp(log_lambda),
    balance_directions, uniform_mixture, counters);
  for (int iter = 0; iter < max_iterations; ++iter) {
    double attempted_proposal_ess = 1.0 / arma::accu(arma::square(current_proposal));
    double attempted_proposal_max = current_proposal.max();
    arma::uvec attempted_positive = arma::find(current_proposal > 0);
    double attempted_proposal_entropy = -arma::accu(
      current_proposal(attempted_positive) % arma::log(current_proposal(attempted_positive)));
    int index = Rcpp::sample(static_cast<int>(gamma.n_elem), 1, false,
                             Rcpp::wrap(current_proposal))[0] - 1;
    bool birth = gamma(index) == 0u;
    arma::uvec proposed = gamma;
    proposed(index) = 1u - proposed(index);
    double proposed_target = exact_log_target(proposed, X, y, a, d, tau, a_phi, b_phi, counters);
    arma::vec proposed_proposal = fast_similarity_probabilities(
      proposed, score_type, gram, cross_y, yty, rss_null, y.n_elem,
      std::exp(log_lambda), balance_directions, uniform_mixture, counters);
    double log_alpha = proposed_target - current_target
      + std::log(proposed_proposal(index)) - std::log(current_proposal(index));
    bool accept = std::log(R::runif(0.0, 1.0)) < std::min(0.0, log_alpha);
    if (accept) {
      gamma = proposed;
      current_target = proposed_target;
      current_proposal = proposed_proposal; // reverse calculation is next forward calculation
      ++accepted; ++window_accepted;
    }
    ++window_attempted;
    window_proposal_ess_sum += attempted_proposal_ess;
    iterations = iter + 1;

    bool lambda_changed = false;
    bool adaptation_active = adaptation_metric == "diversity" || iterations <= burn;
    if (adapt_lambda && adaptation_active && window_attempted == adaptation_window) {
      double rate = static_cast<double>(window_accepted) / window_attempted;
      double step;
      if (adaptation_metric == "diversity") {
        step = adaptation_gain / std::pow(adaptation_offset + updates, adaptation_decay);
        double mean_effective = window_proposal_ess_sum / window_attempted;
        log_lambda += step * (std::log(mean_effective) - std::log(target_effective_neighbors));
      } else {
        step = 0.8 / std::pow(static_cast<double>(updates + 1), 0.6);
        log_lambda += step * (target_acceptance - rate);
      }
      log_lambda = std::max(std::log(lambda_min), std::min(log_lambda, std::log(lambda_max)));
      ++updates; window_accepted = 0; window_attempted = 0;
      window_proposal_ess_sum = 0.0; lambda_changed = true;
    }
    if (lambda_changed) {
      current_proposal = fast_similarity_probabilities(
        gamma, score_type, gram, cross_y, yty, rss_null, y.n_elem,
        std::exp(log_lambda), balance_directions, uniform_mixture, counters);
    }
    if (iterations > burn) {
      if (birth) ++birth_attempts; else ++death_attempts;
      if (accept) {
        ++retained_accepted;
        if (birth) ++birth_accepted; else ++death_accepted;
      }
      proposal_ess_sum += attempted_proposal_ess;
      proposal_entropy_sum += attempted_proposal_entropy;
      proposal_max_sum += attempted_proposal_max;
      pip_sum += arma::conv_to<arma::vec>::from(gamma);
      model_size_trace.push_back(arma::accu(gamma));
      log_target_trace.push_back(current_target);
      lambda_trace.push_back(std::exp(log_lambda));
      if (gamma.n_elem <= 52) {
        double model_index = 0.0;
        for (arma::uword j = 0; j < gamma.n_elem; ++j)
          if (gamma(j)) model_index += std::ldexp(1.0, static_cast<int>(j));
        model_index_trace.push_back(model_index);
      }
    }
    if (iterations >= burn + 50 && elapsed_seconds(total_start) >= time_budget_seconds) break;
  }
  double total_seconds = elapsed_seconds(total_start);
  int retained = model_size_trace.size();
  Rcpp::List output = List::create(
    _["last_gamma"] = gamma,
    _["pips"] = retained ? pip_sum / retained : pip_sum,
    _["model_size"] = Rcpp::wrap(model_size_trace),
    _["log_target"] = Rcpp::wrap(log_target_trace),
    _["lambda_trace"] = Rcpp::wrap(lambda_trace),
    _["final_lambda"] = std::exp(log_lambda),
    _["adaptation_updates"] = updates,
    _["iterations"] = iterations,
    _["retained"] = retained,
    _["acceptance_rate"] = static_cast<double>(accepted) / iterations,
    _["postburn_acceptance_rate"] = retained ? static_cast<double>(retained_accepted) / retained : NA_REAL,
    _["birth_acceptance_rate"] = birth_attempts ? static_cast<double>(birth_accepted) / birth_attempts : NA_REAL,
    _["death_acceptance_rate"] = death_attempts ? static_cast<double>(death_accepted) / death_attempts : NA_REAL,
    _["birth_attempt_fraction"] = retained ? static_cast<double>(birth_attempts) / retained : NA_REAL,
    _["mean_proposal_effective_neighbors"] = retained ? proposal_ess_sum / retained : NA_REAL,
    _["mean_proposal_entropy"] = retained ? proposal_entropy_sum / retained : NA_REAL,
    _["total_seconds"] = total_seconds,
    _["iterations_per_second"] = iterations / total_seconds,
    _["target_seconds"] = counters.target_seconds,
    _["auxiliary_seconds"] = counters.auxiliary_seconds
  );
  output["target_evaluations"] = static_cast<double>(counters.target_evaluations);
  output["neighborhood_factorizations"] = static_cast<double>(counters.ols_fits);
  output["proposal_normalizations"] = static_cast<double>(counters.proposal_normalizations);
  output["model_index"] = Rcpp::wrap(model_index_trace);
  output["mean_max_proposal_probability"] = retained ? proposal_max_sum / retained : NA_REAL;
  output["balanced_directions"] = balance_directions;
  output["adaptation_metric"] = adaptation_metric;
  output["target_effective_neighbors"] = target_effective_neighbors;
  output["uniform_mixture"] = uniform_mixture;
  return output;
}

// [[Rcpp::export]]
Rcpp::List audit_fast_similarity_probabilities(const arma::vec& y,
                                                const arma::mat& predictors,
                                                arma::uvec gamma,
                                                std::string score_type,
                                                double lambda = 0.7,
                                                bool balance_directions = false,
                                                double uniform_mixture = 0.0) {
  arma::mat X = arma::join_rows(arma::ones(y.n_elem), predictors);
  arma::mat gram = X.t() * X;
  arma::vec cross_y = X.t() * y;
  double yty = arma::dot(y, y);
  double rss_null = arma::accu(arma::square(y - arma::mean(y)));
  Counters slow_counters, fast_counters;
  std::string old_method = score_type == "f" ? "f_similarity" : "bic_similarity";
  arma::vec slow = proposal_probabilities(gamma, old_method, X, y, rss_null, lambda,
                                           0.0, 2.0, 2.0, 1.0, 1.0, 3.0,
                                           slow_counters);
  arma::vec fast = fast_similarity_probabilities(gamma, score_type, gram, cross_y,
                                                  yty, rss_null, y.n_elem, lambda,
                                                  balance_directions, uniform_mixture,
                                                  fast_counters);
  return List::create(
    _["max_absolute_difference"] = arma::max(arma::abs(slow - fast)),
    _["slow"] = slow,
    _["fast"] = fast,
    _["slow_candidate_fits"] = static_cast<double>(slow_counters.ols_fits),
    _["fast_factorizations"] = static_cast<double>(fast_counters.ols_fits)
  );
}

// [[Rcpp::export]]
arma::vec evaluate_fast_similarity_probabilities(const arma::vec& y,
                                                  const arma::mat& predictors,
                                                  arma::uvec gamma,
                                                  std::string score_type,
                                                  double lambda = 1.0,
                                                  bool balance_directions = false,
                                                  double uniform_mixture = 0.0) {
  if (score_type != "f" && score_type != "bic" && score_type != "marginal")
    stop("score_type must be f, bic, or marginal");
  arma::mat X = arma::join_rows(arma::ones(y.n_elem), predictors);
  arma::mat gram = X.t()*X;
  arma::vec cross_y = X.t()*y;
  Counters counters;
  return fast_similarity_probabilities(gamma,score_type,gram,cross_y,
    arma::dot(y,y),arma::accu(arma::square(y-arma::mean(y))),y.n_elem,
    lambda,balance_directions,uniform_mixture,counters);
}

// [[Rcpp::export]]
Rcpp::List audit_fast_comparator_probabilities(const arma::vec& y,
                                                const arma::mat& predictors,
                                                arma::uvec gamma,
                                                std::string method,
                                                double a = 2.0, double d = 2.0,
                                                double tau = 1.0,
                                                double a_phi = 1.0,
                                                double b_phi = 3.0) {
  arma::mat X = arma::join_rows(arma::ones(y.n_elem), predictors);
  arma::mat gram = X.t() * X;
  arma::vec cross_y = X.t() * y;
  double yty = arma::dot(y, y);
  double rss_null = arma::accu(arma::square(y - arma::mean(y)));
  Counters slow_counters, fast_counters;
  double current_target = exact_log_target_raw(gamma, X, y, a, d, tau, a_phi, b_phi);
  arma::vec slow = proposal_probabilities(gamma, method, X, y, rss_null, 0.7,
                                           current_target, a, d, tau, a_phi, b_phi,
                                           slow_counters);
  arma::vec fast = fast_comparator_probabilities(gamma, method, gram, cross_y,
                                                  yty, y.n_elem, a, d, tau,
                                                  a_phi, b_phi, fast_counters);
  return List::create(
    _["max_absolute_difference"] = arma::max(arma::abs(slow - fast)),
    _["slow"] = slow,
    _["fast"] = fast,
    _["slow_target_evaluations"] = static_cast<double>(slow_counters.target_evaluations),
    _["fast_factorizations"] = static_cast<double>(fast_counters.ols_fits)
  );
}

// [[Rcpp::export]]
double evaluate_exact_gaussian_log_target(const arma::vec& y,
                                          const arma::mat& predictors,
                                          arma::uvec gamma,
                                          double a = 2.0, double d = 2.0,
                                          double tau = 1.0,
                                          double a_phi = 1.0, double b_phi = 3.0) {
  arma::mat X = arma::join_rows(arma::ones(y.n_elem), predictors);
  return exact_log_target_raw(gamma, X, y, a, d, tau, a_phi, b_phi);
}

// [[Rcpp::export]]
Rcpp::List run_online_gaussian_scaling(const arma::vec& y,
                                       const arma::mat& predictors,
                                       std::string method,
                                       arma::uvec gamma,
                                       double time_budget_seconds,
                                       int burn = 200,
                                       int max_iterations = 1000000,
                                       double lambda = 0.7,
                                       double a = 2.0, double d = 2.0,
                                       double tau = 1.0,
                                       double a_phi = 1.0, double b_phi = 3.0,
                                       std::string proposal_engine = "legacy") {
  RNGScope scope;
  if (time_budget_seconds <= 0 || burn < 0 || max_iterations <= burn)
    stop("Invalid time budget, burn, or maximum iteration count");
  if (predictors.n_rows != y.n_elem || predictors.n_cols != gamma.n_elem)
    stop("Incompatible y, predictors, and gamma dimensions");
  if (proposal_engine != "legacy" && proposal_engine != "rank_one")
    stop("proposal_engine must be legacy or rank_one");

  auto total_start = std::chrono::steady_clock::now();
  arma::mat X = arma::join_rows(arma::ones(y.n_elem), predictors);
  arma::mat gram = X.t() * X;
  arma::vec cross_y = X.t() * y;
  double yty = arma::dot(y, y);
  double rss_null = arma::accu(arma::square(y - arma::mean(y)));
  Counters counters;
  arma::vec pip_sum(gamma.n_elem, arma::fill::zeros);
  std::vector<double> model_size_trace;
  std::vector<double> log_target_trace;
  model_size_trace.reserve(std::min(max_iterations - burn, 200000));
  log_target_trace.reserve(std::min(max_iterations - burn, 200000));
  int accepted = 0;
  int iterations = 0;
  double current_target = exact_log_target(gamma, X, y, a, d, tau, a_phi, b_phi, counters);
  auto forward_start = std::chrono::steady_clock::now();
  arma::vec current_proposal = proposal_engine == "rank_one"
    ? fast_comparator_probabilities(gamma, method, gram, cross_y, yty, y.n_elem,
                                    a, d, tau, a_phi, b_phi, counters)
    : proposal_probabilities(gamma, method, X, y, rss_null, lambda, current_target,
                             a, d, tau, a_phi, b_phi, counters);
  counters.forward_seconds += elapsed_seconds(forward_start);

  for (int iter = 0; iter < max_iterations; ++iter) {
    Rcpp::NumericVector proposal_weights = Rcpp::wrap(current_proposal);
    int index = Rcpp::sample(static_cast<int>(gamma.n_elem), 1, false,
                             proposal_weights)[0] - 1;
    arma::uvec proposed = gamma;
    proposed(index) = 1u - proposed(index);
    double proposed_target = exact_log_target(proposed, X, y, a, d, tau, a_phi, b_phi, counters);

    auto reverse_start = std::chrono::steady_clock::now();
    arma::vec proposed_proposal = proposal_engine == "rank_one"
      ? fast_comparator_probabilities(proposed, method, gram, cross_y, yty, y.n_elem,
                                      a, d, tau, a_phi, b_phi, counters)
      : proposal_probabilities(proposed, method, X, y, rss_null, lambda, proposed_target,
                               a, d, tau, a_phi, b_phi, counters);
    counters.reverse_seconds += elapsed_seconds(reverse_start);
    double log_alpha = proposed_target - current_target
      + std::log(proposed_proposal(index)) - std::log(current_proposal(index));
    if (std::log(R::runif(0.0, 1.0)) < std::min(0.0, log_alpha)) {
      gamma = proposed;
      current_target = proposed_target;
      current_proposal = proposed_proposal;
      accepted += 1;
    }
    iterations = iter + 1;
    if (iterations > burn) {
      pip_sum += arma::conv_to<arma::vec>::from(gamma);
      model_size_trace.push_back(arma::accu(gamma));
      log_target_trace.push_back(current_target);
    }
    if (iterations >= burn + 50 && elapsed_seconds(total_start) >= time_budget_seconds)
      break;
  }

  double total_seconds = elapsed_seconds(total_start);
  int retained = model_size_trace.size();
  arma::vec pips = retained > 0 ? pip_sum / retained : pip_sum;
  return List::create(
    _["last_gamma"] = gamma,
    _["pips"] = pips,
    _["model_size"] = Rcpp::wrap(model_size_trace),
    _["log_target"] = Rcpp::wrap(log_target_trace),
    _["iterations"] = iterations,
    _["retained"] = retained,
    _["hit_iteration_cap"] = iterations >= max_iterations,
    _["acceptance_rate"] = static_cast<double>(accepted) / iterations,
    _["total_seconds"] = total_seconds,
    _["iterations_per_second"] = iterations / total_seconds,
    _["target_seconds"] = counters.target_seconds,
    _["auxiliary_seconds"] = counters.auxiliary_seconds,
    _["forward_seconds"] = counters.forward_seconds,
    _["reverse_seconds"] = counters.reverse_seconds,
    _["target_evaluations"] = static_cast<double>(counters.target_evaluations),
    _["ols_fits"] = static_cast<double>(counters.ols_fits),
    _["proposal_normalizations"] = static_cast<double>(counters.proposal_normalizations)
  );
}

// [[Rcpp::export]]
Rcpp::List run_online_gaussian_kernel(const arma::vec& y,
                                      const arma::mat& predictors,
                                      std::string method,
                                      arma::uvec gamma,
                                      int runs, int burn,
                                      double lambda = 0.7,
                                      double a = 2.0, double d = 2.0,
                                      double tau = 1.0,
                                      double a_phi = 1.0, double b_phi = 3.0) {
  RNGScope scope;
  if (burn < 0 || runs <= burn) stop("Require 0 <= burn < runs");
  if (predictors.n_rows != y.n_elem || predictors.n_cols != gamma.n_elem)
    stop("Incompatible y, predictors, and gamma dimensions");
  arma::mat X = arma::join_rows(arma::ones(y.n_elem), predictors);
  double rss_null = arma::accu(arma::square(y - arma::mean(y)));
  Counters counters;
  arma::umat draws(runs - burn, gamma.n_elem, arma::fill::zeros);
  arma::vec model_size(runs - burn, arma::fill::zeros);
  arma::vec log_target_trace(runs - burn, arma::fill::zeros);
  int accepted = 0;
  double current_target = exact_log_target(gamma, X, y, a, d, tau, a_phi, b_phi, counters);
  auto total_start = std::chrono::steady_clock::now();

  for (int iter = 0; iter < runs; ++iter) {
    auto forward_start = std::chrono::steady_clock::now();
    arma::vec forward = proposal_probabilities(gamma, method, X, y, rss_null, lambda,
                                                current_target, a, d, tau, a_phi, b_phi, counters);
    counters.forward_seconds += elapsed_seconds(forward_start);
    Rcpp::NumericVector proposal_weights = Rcpp::wrap(forward);
    int index = Rcpp::sample(static_cast<int>(gamma.n_elem), 1, false,
                             proposal_weights)[0] - 1;
    arma::uvec proposed = gamma;
    proposed(index) = 1u - proposed(index);
    double proposed_target = exact_log_target(proposed, X, y, a, d, tau, a_phi, b_phi, counters);

    auto reverse_start = std::chrono::steady_clock::now();
    arma::vec reverse = proposal_probabilities(proposed, method, X, y, rss_null, lambda,
                                                proposed_target, a, d, tau, a_phi, b_phi, counters);
    counters.reverse_seconds += elapsed_seconds(reverse_start);
    double log_alpha = proposed_target - current_target
      + std::log(reverse(index)) - std::log(forward(index));
    if (std::log(R::runif(0.0, 1.0)) < std::min(0.0, log_alpha)) {
      gamma = proposed;
      current_target = proposed_target;
      accepted += 1;
    }
    if (iter >= burn) {
      int out = iter - burn;
      draws.row(out) = gamma.t();
      model_size(out) = arma::accu(gamma);
      log_target_trace(out) = current_target;
    }
  }
  double total_seconds = elapsed_seconds(total_start);
  return List::create(
    _["gamma"] = draws,
    _["model_size"] = model_size,
    _["log_target"] = log_target_trace,
    _["acceptance_rate"] = static_cast<double>(accepted) / runs,
    _["total_seconds"] = total_seconds,
    _["seconds_per_iteration"] = total_seconds / runs,
    _["target_seconds"] = counters.target_seconds,
    _["auxiliary_seconds"] = counters.auxiliary_seconds,
    _["forward_seconds"] = counters.forward_seconds,
    _["reverse_seconds"] = counters.reverse_seconds,
    _["target_evaluations"] = static_cast<double>(counters.target_evaluations),
    _["ols_fits"] = static_cast<double>(counters.ols_fits),
    _["proposal_normalizations"] = static_cast<double>(counters.proposal_normalizations)
  );
}

// Fixed-length comparator used by the production exact-posterior audit.  The
// reverse proposal evaluated at an accepted state is cached as the next
// forward proposal, so timing includes neighborhood scoring without doing
// avoidable duplicate work.
// [[Rcpp::export]]
Rcpp::List run_exact_comparator_kernel(const arma::vec& y,
                                       const arma::mat& predictors,
                                       std::string method,
                                       arma::uvec gamma,
                                       int runs, int burn,
                                       double a = 2.0, double d = 2.0,
                                       double tau = 1.0,
                                       double a_phi = 1.0, double b_phi = 3.0,
                                       std::string proposal_engine = "legacy") {
  RNGScope scope;
  if (burn < 0 || runs <= burn) stop("Require 0 <= burn < runs");
  if (method != "uniform" && method != "exact_locally_balanced" &&
      method != "bic_locally_balanced")
    stop("Comparator method must be uniform, exact_locally_balanced, or bic_locally_balanced");
  if (proposal_engine != "legacy" && proposal_engine != "rank_one")
    stop("proposal_engine must be legacy or rank_one");
  if (predictors.n_rows != y.n_elem || predictors.n_cols != gamma.n_elem)
    stop("Incompatible y, predictors, and gamma dimensions");

  arma::mat X = arma::join_rows(arma::ones(y.n_elem), predictors);
  arma::mat gram = X.t() * X;
  arma::vec cross_y = X.t() * y;
  double yty = arma::dot(y, y);
  double rss_null = arma::accu(arma::square(y - arma::mean(y)));
  Counters counters;
  arma::umat draws(runs - burn, gamma.n_elem, arma::fill::zeros);
  arma::vec model_size(runs - burn, arma::fill::zeros);
  arma::vec log_target_trace(runs - burn, arma::fill::zeros);
  int accepted = 0;
  auto total_start = std::chrono::steady_clock::now();
  double current_target = exact_log_target(gamma, X, y, a, d, tau,
                                            a_phi, b_phi, counters);
  arma::vec current_proposal = proposal_engine == "rank_one"
    ? fast_comparator_probabilities(gamma, method, gram, cross_y, yty, y.n_elem,
                                    a, d, tau, a_phi, b_phi, counters)
    : proposal_probabilities(gamma, method, X, y, rss_null, 0.7, current_target,
                             a, d, tau, a_phi, b_phi, counters);
  for (int iter = 0; iter < runs; ++iter) {
    Rcpp::NumericVector proposal_weights = Rcpp::wrap(current_proposal);
    int index = Rcpp::sample(static_cast<int>(gamma.n_elem), 1, false,
                             proposal_weights)[0] - 1;
    arma::uvec proposed = gamma;
    proposed(index) = 1u - proposed(index);
    double proposed_target = exact_log_target(proposed, X, y, a, d, tau,
                                               a_phi, b_phi, counters);

    auto reverse_start = std::chrono::steady_clock::now();
    arma::vec proposed_proposal = proposal_engine == "rank_one"
      ? fast_comparator_probabilities(proposed, method, gram, cross_y, yty, y.n_elem,
                                      a, d, tau, a_phi, b_phi, counters)
      : proposal_probabilities(proposed, method, X, y, rss_null, 0.7, proposed_target,
                               a, d, tau, a_phi, b_phi, counters);
    counters.reverse_seconds += elapsed_seconds(reverse_start);
    double log_alpha = proposed_target - current_target
      + std::log(proposed_proposal(index)) - std::log(current_proposal(index));
    if (std::log(R::runif(0.0, 1.0)) < std::min(0.0, log_alpha)) {
      gamma = proposed;
      current_target = proposed_target;
      current_proposal = proposed_proposal;
      accepted += 1;
    }
    if (iter >= burn) {
      int out = iter - burn;
      draws.row(out) = gamma.t();
      model_size(out) = arma::accu(gamma);
      log_target_trace(out) = current_target;
    }
  }
  double total_seconds = elapsed_seconds(total_start);
  return List::create(
    _["gamma"] = draws,
    _["model_size"] = model_size,
    _["log_target"] = log_target_trace,
    _["acceptance_rate"] = static_cast<double>(accepted) / runs,
    _["total_seconds"] = total_seconds,
    _["seconds_per_iteration"] = total_seconds / runs,
    _["target_seconds"] = counters.target_seconds,
    _["auxiliary_seconds"] = counters.auxiliary_seconds,
    _["reverse_seconds"] = counters.reverse_seconds,
    _["target_evaluations"] = static_cast<double>(counters.target_evaluations),
    _["ols_fits"] = static_cast<double>(counters.ols_fits),
    _["proposal_normalizations"] = static_cast<double>(counters.proposal_normalizations)
  );
}
