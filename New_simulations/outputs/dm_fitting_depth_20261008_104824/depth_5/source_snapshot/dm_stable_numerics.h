#ifndef DM_STABLE_NUMERICS_H
#define DM_STABLE_NUMERICS_H

#include <cmath>
#include <limits>
#include <atomic>

namespace dm_numeric_diagnostics {
  static std::atomic<int> rejected_ratios{0};
  static std::atomic<int> precision_projections{0};
  static std::atomic<int> optimizer_fallbacks{0};
}

// Integer-count rising factorials evaluated on the log-concentration scale.
// No clipping of predictors or changes to the statistical likelihood.
struct DMRising {
  double correction; // log((a)_m) - m*log(a)
  double first;      // a * [digamma(a+m)-digamma(a)]
  double second;     // derivative of first with respect to log(a)
};

inline DMRising dm_rising(double log_a, int count) {
  if (count == 0) return {0.0, 0.0, 0.0};
  if (!std::isfinite(log_a) || count < 0) {
    double bad = std::numeric_limits<double>::quiet_NaN();
    return {bad, bad, bad};
  }
  const double m = count;
  if (log_a < -36.0) {
    return {(1.0-m)*log_a + R::lgammafn(m), 1.0, 0.0};
  }
  // Convergent expansion when every k/a is small.  Eight terms bound the
  // neglected relative correction below roundoff at this threshold.
  if (log_a > std::log(std::max(1.0, m)) + std::log(100.0)) {
    const long double inv_a = std::exp(-static_cast<long double>(log_a));
    long double powers[9] = {m, 0, 0, 0, 0, 0, 0, 0, 0};
    // Faulhaber's recurrence: sum_{k=0}^{m-1} k^r.
    for (int r=1; r<=8; ++r) {
      long double sum=0, choose=1;
      for (int j=0; j<r; ++j) {
        sum += choose*powers[j];
        choose *= static_cast<long double>(r+1-j)/(j+1);
      }
      powers[r]=(std::pow(static_cast<long double>(m),r+1)-sum)/(r+1);
    }
    long double correction=0, first=m, second=0, scale=1;
    for (int r=1; r<=8; ++r) {
      scale *= inv_a;
      const long double term=powers[r]*scale;
      const int sign=(r%2) ? 1 : -1;
      correction += sign*term/r;
      first -= sign*term;
      second += sign*r*term;
    }
    return {static_cast<double>(correction),static_cast<double>(first),
            static_cast<double>(second)};
  }
  const double a=std::exp(log_a);
  if (count <= 64) {
    double correction=0, first=0, second=0;
    for (int k=0; k<count; ++k) {
      const double ratio=a/(a+k);
      correction += std::log1p(k/a);
      first += ratio;
      second += ratio*(1-ratio);
    }
    return {correction,first,second};
  }
  const double first=a*(R::digamma(a+m)-R::digamma(a));
  const double second=first+a*a*(R::trigamma(a+m)-R::trigamma(a));
  return {R::lgammafn(m)-R::lbeta(a,m)-m*log_a,first,second};
}

inline double dm_log_sum_exp(const arma::rowvec& eta) {
  if (!eta.is_finite()) return std::numeric_limits<double>::quiet_NaN();
  const double top=eta.max();
  return top+std::log(arma::accu(arma::exp(eta-top)));
}

inline double dm_row_log_kernel(const arma::irowvec& y,
                                const arma::rowvec& eta) {
  const double log_sum=dm_log_sum_exp(eta);
  if (!std::isfinite(log_sum)) return -std::numeric_limits<double>::infinity();
  const int total=arma::accu(y);
  long double value=-dm_rising(log_sum,total).correction;
  for (arma::uword j=0;j<y.n_elem;++j) {
    if (y(j)>0) {
      value += static_cast<long double>(y(j))*(eta(j)-log_sum)
        +dm_rising(eta(j),y(j)).correction;
    }
  }
  return static_cast<double>(value);
}

inline double dm_stable_log_kernel(const arma::imat& y,const arma::mat& eta) {
  if (!eta.is_finite()) return -std::numeric_limits<double>::infinity();
  long double value=0;
  for (arma::uword i=0;i<y.n_rows;++i) {
    const double row=dm_row_log_kernel(y.row(i),eta.row(i));
    if (!std::isfinite(row)) return -std::numeric_limits<double>::infinity();
    value += row;
  }
  return static_cast<double>(value);
}

inline double dm_stable_replacement_delta(const arma::imat& y,
    const arma::mat& alpha_mat,const arma::mat& xbeta,int j,
    const arma::vec& replacement) {
  if (!replacement.is_finite()) return -std::numeric_limits<double>::infinity();
  long double delta=0;
  for (arma::uword i=0;i<y.n_rows;++i) {
    arma::rowvec old_eta=alpha_mat.row(i)+xbeta.row(i);
    arma::rowvec new_eta=old_eta;
    new_eta(j)=alpha_mat(i,j)+replacement(i);
    const double old_sum=dm_log_sum_exp(old_eta);
    const double new_sum=dm_log_sum_exp(new_eta);
    if (!std::isfinite(old_sum)||!std::isfinite(new_sum))
      return -std::numeric_limits<double>::infinity();
    const int total=arma::accu(y.row(i));
    delta += static_cast<long double>(y(i,j))*(new_eta(j)-old_eta(j))
      -static_cast<long double>(total)*(new_sum-old_sum)
      +dm_rising(new_eta(j),y(i,j)).correction
      -dm_rising(old_eta(j),y(i,j)).correction
      -dm_rising(new_sum,total).correction
      +dm_rising(old_sum,total).correction;
  }
  return static_cast<double>(delta);
}

inline double dm_safe_log_acceptance(double ratio) {
  // Reject nonfinite ratios, including +Inf: a finite current state and
  // finite proposal components must never produce an infinite MH ratio.
  if (std::isfinite(ratio)) return std::min(0.0,ratio);
  ++dm_numeric_diagnostics::rejected_ratios;
  return -std::numeric_limits<double>::infinity();
}

inline arma::mat dm_gaussian_precision(const arma::mat& raw) {
  if (!raw.is_finite()) Rcpp::stop("Nonfinite DM coefficient precision");
  arma::mat q=arma::symmatu(raw), chol;
  if (arma::chol(chol,q)) return q;
  ++dm_numeric_diagnostics::precision_projections;
  // Finite indefinite curvature can arise at an approximate mode.  Project
  // the proposal precision onto the SPD cone, avoiding the old 1e-8 floor
  // (which allowed standard deviations of 10000).  This changes only q.
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  if (!arma::eig_sym(eigenvalues,eigenvectors,q))
    Rcpp::stop("DM precision eigendecomposition failed");
  eigenvalues=arma::clamp(eigenvalues,0.1,std::numeric_limits<double>::infinity());
  q=arma::symmatu(eigenvectors*arma::diagmat(eigenvalues)*eigenvectors.t());
  if (!arma::chol(chol,q)) Rcpp::stop("DM stabilized precision is not SPD");
  return q;
}
#endif
