#ifndef GAUSSIAN_REVISION_SHARED_H
#define GAUSSIAN_REVISION_SHARED_H
#include <limits>

// Shared numerical conventions for the production and benchmark kernels.
namespace gaussrev {
inline double prior(double k, double p, double a, double b) {
  return R::lbeta(a + k, b + p - k) - R::lbeta(a, b);
}
inline double rss(const arma::mat& gram, const arma::vec& rhs, double yty) {
  arma::mat inverse;
  arma::vec beta;
  if (arma::rcond(gram) > 1e-12 && arma::inv_sympd(inverse, gram)) beta = inverse * rhs;
  else beta = arma::pinv(gram) * rhs;
  return std::max(yty - arma::dot(rhs, beta), 1e-14);
}
inline double tail_score(double logp) {
  return std::max(-std::max(logp, std::log(std::numeric_limits<double>::min())) /
                  std::log(10.0), 0.0);
}
inline double f_score(int k, int n, double rss0, double residual) {
  // A neutral finite auxiliary score when residual degrees of freedom do not
  // exist; this does NOT restrict the conjugate posterior's model support.
  if (k == 0 || n - k - 1 <= 0) return 0.0;
  double statistic = std::max(((rss0 - residual) / k) /
                              (residual / (n - k - 1)), 0.0);
  return tail_score(R::pf(statistic, 1, n - k - 1, false, true));
}
inline double lr_score(int k, int n, double rss0, double residual,
                       const std::string& type) {
  if (k == 0) return 0.0;
  double statistic = n * (std::log(std::max(rss0, 1e-14)) - std::log(residual));
  if (type == "bic") statistic = std::max(statistic - std::log((double)n) * k, 0.0);
  else statistic = std::max(statistic / std::sqrt((double)k), 0.0);
  return tail_score(R::pchisq(statistic, 1, false, true));
}
}
#endif
