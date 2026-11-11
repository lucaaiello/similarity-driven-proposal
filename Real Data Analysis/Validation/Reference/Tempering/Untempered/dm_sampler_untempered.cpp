// [[Rcpp::plugins(cpp14,openmp)]]
#include <algorithm>
#include <chrono>
#include <cstdio>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif
#include <iomanip>
#include <RcppEnsmallen.h>
#include <Rcpp.h>
#include <omp.h>
#include <sstream>
#include <unordered_map>
#include <memory>
#include <cstdint>
#include "dm_stable_numerics.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppEnsmallen)]]
// [[Rcpp::depends(RcppArmadillo)]]

////////////////////////
// Optimization stuff //
////////////////////////

// V2 computational copy. Production, archived baseline and tested V1 are immutable.
// All statistical kernels/settings stay intact; switches below are diagnostics only.
namespace dm_runtime {
  enum Metric { construction=0, objective=1, hessian=2, fitting=3, neighborhood=4,
    replacement_context=5, replacement=6, context_key=7, metric_count=8 };
  static bool enabled=false;
  static bool use_cached=true;
  static bool cache_replacements=true;
  static bool split_numerics=true;
  static bool share_neighborhood_context=true;
  static std::uint64_t generation=0;
  static std::atomic<std::uint64_t> calls[metric_count]{};
  static std::atomic<std::uint64_t> nanoseconds[metric_count]{};
  static std::atomic<std::uint64_t> fit_cache_hits{0};
  static std::atomic<std::uint64_t> conditional_context_hits{0};
  static std::atomic<std::uint64_t> conditional_context_misses{0};
  struct Scope {
    Metric metric;
    bool active;
    std::chrono::steady_clock::time_point start;
    explicit Scope(Metric metric):metric(metric),active(enabled) {
      if(active) { calls[metric].fetch_add(1,std::memory_order_relaxed); start=std::chrono::steady_clock::now(); }
    }
    ~Scope() {
      if(active) nanoseconds[metric].fetch_add(static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-start).count()),
        std::memory_order_relaxed);
    }
  };
}

// [[Rcpp::export]]
void dm_runtime_v2_options(bool cache_replacements=true,bool split_numerics=true,
    bool share_neighborhood_context=true) {
  dm_runtime::cache_replacements=cache_replacements;
  dm_runtime::split_numerics=split_numerics;
  dm_runtime::share_neighborhood_context=share_neighborhood_context;
  ++dm_runtime::generation;
}

// [[Rcpp::export]]
void dm_runtime_use_cached(bool use_cached=true) {
  dm_runtime::use_cached=use_cached;
  // Engine switches invalidate caches, even if the requested engine is unchanged.
  ++dm_runtime::generation;
}

// [[Rcpp::export]]
void dm_runtime_profile_reset(bool enabled=true) {
  dm_runtime::enabled=enabled;
  for(int k=0;k<dm_runtime::metric_count;++k) {dm_runtime::calls[k]=0;dm_runtime::nanoseconds[k]=0;}
  dm_runtime::fit_cache_hits=0;
  dm_runtime::conditional_context_hits=0;
  dm_runtime::conditional_context_misses=0;
}

// [[Rcpp::export]]
Rcpp::List dm_runtime_profile_snapshot() {
  Rcpp::CharacterVector labels=Rcpp::CharacterVector::create("constructor","objective","hessian","fit","neighborhood",
    "replacement_context","replacement","context_key");
  Rcpp::NumericVector count(dm_runtime::metric_count),ns(dm_runtime::metric_count);
  for(int k=0;k<dm_runtime::metric_count;++k) { count[k]=dm_runtime::calls[k].load();ns[k]=dm_runtime::nanoseconds[k].load(); }
  return Rcpp::List::create(Rcpp::Named("metrics")=Rcpp::DataFrame::create(
      Rcpp::Named("metric")=labels,Rcpp::Named("calls")=count,Rcpp::Named("total_nanoseconds")=ns),
    Rcpp::Named("enabled")=dm_runtime::enabled,
    Rcpp::Named("engine")=dm_runtime::use_cached ? "conditional_cached_v2" : "full_reference",
    Rcpp::Named("cache_replacements")=dm_runtime::cache_replacements,
    Rcpp::Named("split_numerics")=dm_runtime::split_numerics,
    Rcpp::Named("share_neighborhood_context")=dm_runtime::share_neighborhood_context,
    Rcpp::Named("fit_cache_hits")=static_cast<double>(dm_runtime::fit_cache_hits.load()),
    Rcpp::Named("conditional_context_hits")=static_cast<double>(dm_runtime::conditional_context_hits.load()),
    Rcpp::Named("conditional_context_misses")=static_cast<double>(dm_runtime::conditional_context_misses.load()),
    Rcpp::Named("timing_note")="Nested, potentially overlapping summed work times; not additive wall-clock components.");
}

class DMRegressionFunctionFullReference {
public:
  DMRegressionFunctionFullReference(const arma::mat& X_j, 
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
    
    if (j_aux < 0 || j_aux >= J || intercepts.n_elem != Y.n_cols ||
        Xbeta.n_rows != Y.n_rows || Xbeta.n_cols != Y.n_cols ||
        X_j.n_rows != Y.n_rows || PP.n_elem != Y.n_cols)
      Rcpp::stop("Invalid dimensions in DM conditional objective");
    P_j_aux = PP(j_aux);
    if (P_j_aux != static_cast<int>(X_j.n_cols))
      Rcpp::stop("DM conditional design and coefficient count disagree");
    
  }
  
  //////////////////////////////////////
  // Objective and Gradient
  //////////////////////////////////////
  
  double EvaluateWithGradient(const arma::mat& beta_j, arma::mat& grad_flat) {
    grad_flat.zeros(P_j_aux,1);
    if (P_j_aux==0) return 0.0;
    if (!beta_j.is_finite()) return std::numeric_limits<double>::infinity();
    const int threads=omp_get_max_threads();
    std::vector<arma::vec> gradients(threads,arma::vec(P_j_aux,arma::fill::zeros));
    std::vector<double> values(threads,0.0);
#pragma omp parallel for schedule(static)
    for (int i=0;i<static_cast<int>(n);++i) {
      int tid=omp_get_thread_num();
      arma::rowvec eta=intercepts+Xbeta.row(i);
      eta(j_aux)=intercepts(j_aux)+arma::dot(X_j.row(i),beta_j.col(0));
      double log_sum=dm_log_sum_exp(eta);
      double ll=dm_row_log_kernel(Y.row(i),eta);
      if (!std::isfinite(ll)) {
        values[tid]=std::numeric_limits<double>::infinity();
        continue;
      }
      DMRising own=dm_rising(eta(j_aux),Y(i,j_aux));
      DMRising total=dm_rising(log_sum,arma::accu(Y.row(i)));
      double weight=std::exp(eta(j_aux)-log_sum);
      gradients[tid] -= (own.first-weight*total.first)*X_j.row(i).t();
      values[tid] -= ll;
    }
    double value=0;
    for (int t=0;t<threads;++t) { grad_flat+=gradients[t]; value+=values[t]; }
    grad_flat+=lambda*beta_j;
    value+=0.5*lambda*arma::accu(arma::square(beta_j));
    return value;
  }

  arma::mat ComputeHessian(const arma::mat& beta_j) {
    if (P_j_aux==0) return arma::mat(0,0);
    const int threads=omp_get_max_threads();
    std::vector<arma::mat> blocks(threads,arma::mat(P_j_aux,P_j_aux,arma::fill::zeros));
#pragma omp parallel for schedule(static)
    for (int i=0;i<static_cast<int>(n);++i) {
      int tid=omp_get_thread_num();
      arma::rowvec eta=intercepts+Xbeta.row(i);
      eta(j_aux)=intercepts(j_aux)+arma::dot(X_j.row(i),beta_j.col(0));
      double log_sum=dm_log_sum_exp(eta);
      DMRising own=dm_rising(eta(j_aux),Y(i,j_aux));
      DMRising total=dm_rising(log_sum,arma::accu(Y.row(i)));
      double weight=std::exp(eta(j_aux)-log_sum);
      double curvature=own.second-weight*total.first+
        weight*weight*(total.first-total.second);
      blocks[tid] += curvature*(X_j.row(i).t()*X_j.row(i));
    }
    arma::mat h(P_j_aux,P_j_aux,arma::fill::zeros);
    for (int t=0;t<threads;++t) h+=blocks[t];
    h.diag()-=lambda;
    return arma::symmatu(h);
  }
  
private:
  const arma::mat& X_j;
  const arma::imat& Y;
  // Own this row: callers may pass alpha.t(), whose evaluated temporary
  // expires at the end of the constructor statement in audit entry points.
  const arma::rowvec intercepts;
  const arma::ivec& PP;
  const arma::mat& Xbeta;
  const int j_aux;
  const double lambda;
  size_t n;
  int J;
  int P_j_aux;
};

// Exact-context identity. Current genus Xbeta is deliberately excluded because
// its entire predictor is replaced by intercept_j + X_j*beta in every evaluation.
// Intercepts (INCLUDING j), all other genera, counts and engine generation matter.
std::string dm_runtime_context_key(const arma::imat& Y,const arma::rowvec& alpha,
    const arma::mat& Xbeta,int j) {
  dm_runtime::Scope timer(dm_runtime::context_key);
  std::string context;
  auto append=[&context](const void* data,size_t bytes) {
    if(bytes) context.append(static_cast<const char*>(data),bytes);
  };
  append(&dm_runtime::generation,sizeof(dm_runtime::generation));
  append(&j,sizeof(j));
  const arma::uword shape[]={Y.n_rows,Y.n_cols,alpha.n_elem};
  append(shape,sizeof(shape));
  append(Y.memptr(),Y.n_elem*sizeof(arma::sword));
  append(alpha.memptr(),alpha.n_elem*sizeof(double));
  arma::mat rest=Xbeta;rest.col(j).zeros();
  append(rest.memptr(),rest.n_elem*sizeof(double));
  return context;
}

// Dispatch retains the untouched full helper as an independently selectable path.
inline double dm_runtime_value(double eta,int count) {
  return dm_runtime::use_cached && dm_runtime::split_numerics ?
    dm_rising_value(eta,count) : dm_rising(eta,count).correction;
}
inline DMRising dm_runtime_value_first(double eta,int count) {
  return dm_runtime::use_cached && dm_runtime::split_numerics ?
    dm_rising_value_first(eta,count) : dm_rising(eta,count);
}
inline DMRising dm_runtime_derivatives(double eta,int count) {
  return dm_runtime::use_cached && dm_runtime::split_numerics ?
    dm_rising_derivatives(eta,count) : dm_rising(eta,count);
}

// Immutable, stack-scoped to ONE physical neighbourhood. Unlike fitting context,
// this includes the current own-genus predictor: the delta is relative to it.
// Nothing persists after a conditioning state changes; no history/stale scores.
struct DMRuntimeReplacementContext {
  arma::vec old_own,old_sum,log_other,intercept,old_own_value,old_total_value;
  arma::ivec own_count,total;
  bool valid,single;
  DMRuntimeReplacementContext(const arma::imat& Y,const arma::mat& alpha_mat,
      const arma::mat& xb,int j):old_own(Y.n_rows),old_sum(Y.n_rows),log_other(Y.n_rows),
      intercept(alpha_mat.col(j)),old_own_value(Y.n_rows),old_total_value(Y.n_rows),
      own_count(Y.col(j)),total(Y.n_rows),valid(true),single(Y.n_cols==1) {
    dm_runtime::Scope timer(dm_runtime::replacement_context);
    for(arma::uword i=0;i<Y.n_rows;++i) {
      arma::rowvec eta=alpha_mat.row(i)+xb.row(i);
      total(i)=arma::accu(Y.row(i));old_own(i)=eta(j);
      old_sum(i)=dm_log_sum_exp(eta);
      if(!std::isfinite(old_sum(i)))valid=false;
      if(!single) {
        double top=-std::numeric_limits<double>::infinity();
        for(arma::uword g=0;g<Y.n_cols;++g)if(static_cast<int>(g)!=j)top=std::max(top,eta(g));
        double mass=0;
        for(arma::uword g=0;g<Y.n_cols;++g)if(static_cast<int>(g)!=j)mass+=std::exp(eta(g)-top);
        log_other(i)=top+std::log(mass);
        if(!std::isfinite(log_other(i)))valid=false;
      }
      old_own_value(i)=dm_runtime_value(old_own(i),own_count(i));
      old_total_value(i)=dm_runtime_value(old_sum(i),total(i));
    }
  }
  double delta(const arma::vec& replacement)const {
    dm_runtime::Scope timer(dm_runtime::replacement);
    if(!valid || replacement.n_elem!=old_own.n_elem || !replacement.is_finite())
      return -std::numeric_limits<double>::infinity();
    long double value=0;
    for(arma::uword i=0;i<replacement.n_elem;++i) {
      const double own=intercept(i)+replacement(i);
      double sum=own;
      if(!single) {
        const double lm=log_other(i),top=std::max(own,lm);
        sum=top+std::log(std::exp(own-top)+std::exp(lm-top));
      }
      if(!std::isfinite(sum))return -std::numeric_limits<double>::infinity();
      value+=static_cast<long double>(own_count(i))*(own-old_own(i))
        -static_cast<long double>(total(i))*(sum-old_sum(i))
        +dm_runtime_value(own,own_count(i))-old_own_value(i)
        -dm_runtime_value(sum,total(i))+old_total_value(i);
    }
    return static_cast<double>(value);
  }
};

struct DMCachedConditionalContext {
  arma::vec log_other;
  arma::ivec total,other_count,own_count;
  std::vector<long double> fixed;
  double own_intercept;
  bool single_category,valid;
  DMCachedConditionalContext(const arma::imat& Y,const arma::rowvec& alpha,
      const arma::mat& Xbeta,int j):log_other(Y.n_rows),total(Y.n_rows),
      other_count(Y.n_rows),own_count(Y.col(j)),fixed(Y.n_rows,0),
      own_intercept(alpha(j)),single_category(Y.n_cols==1),valid(std::isfinite(alpha(j))) {
    for(arma::uword i=0;i<Y.n_rows;++i) {
      total(i)=arma::accu(Y.row(i));other_count(i)=total(i)-own_count(i);
      if(single_category) {log_other(i)=0;continue;}
      double top=-std::numeric_limits<double>::infinity();
      for(arma::uword g=0;g<Y.n_cols;++g)if(static_cast<int>(g)!=j) {
        const double value=alpha(g)+Xbeta(i,g);
        if(!std::isfinite(value))valid=false;
        top=std::max(top,value);
      }
      double mass=0;
      for(arma::uword g=0;g<Y.n_cols;++g)if(static_cast<int>(g)!=j)
        mass+=std::exp(alpha(g)+Xbeta(i,g)-top);
      const double lm=top+std::log(mass);log_other(i)=lm;
      if(!std::isfinite(lm))valid=false;
      long double constant=0;
      for(arma::uword g=0;g<Y.n_cols;++g)if(static_cast<int>(g)!=j && Y(i,g)>0) {
        const double value=alpha(g)+Xbeta(i,g);
        constant+=static_cast<long double>(Y(i,g))*(value-lm)+dm_runtime_value(value,Y(i,g));
      }
      fixed[i]=constant;
    }
  }
};

std::shared_ptr<const DMCachedConditionalContext> dm_runtime_conditional_context(
    const arma::imat& Y,const arma::rowvec& alpha,const arma::mat& Xbeta,int j,
    const std::string* supplied_key) {
  struct Cache { std::string key;std::shared_ptr<const DMCachedConditionalContext> data; };
  static thread_local Cache cache;
  const std::string generated=supplied_key ? std::string() : dm_runtime_context_key(Y,alpha,Xbeta,j);
  const std::string& key=supplied_key ? *supplied_key : generated;
  if(cache.data && cache.key==key) {
    if(dm_runtime::enabled)++dm_runtime::conditional_context_hits;
    return cache.data;
  }
  if(dm_runtime::enabled)++dm_runtime::conditional_context_misses;
  auto result=std::make_shared<DMCachedConditionalContext>(Y,alpha,Xbeta,j);
  cache.key=key;cache.data=result;
  return result;
}

class DMRegressionFunctionCached {
public:
  DMRegressionFunctionCached(const arma::mat& design,const arma::imat& Y,
      const arma::rowvec& alpha,const arma::ivec& PP,const arma::mat& Xbeta,int j,
      double penalty,const std::string* key,
      std::shared_ptr<const DMCachedConditionalContext> shared_context=nullptr):
      X_j(design),lambda(penalty),n(Y.n_rows),P(PP(j)),
      context(shared_context ? shared_context : dm_runtime_conditional_context(Y,alpha,Xbeta,j,key)) {}

  double EvaluateWithGradient(const arma::mat& beta,arma::mat& gradient) {
    gradient.zeros(P,1);
    if(P==0)return 0.0; // preserve original empty-design sentinel
    if(!beta.is_finite() || !context->valid)return std::numeric_limits<double>::infinity();
    if(context->single_category) {
      arma::vec eta=X_j*beta;eta+=context->own_intercept;
      if(!eta.is_finite())return std::numeric_limits<double>::infinity();
      gradient=lambda*beta;return .5*lambda*arma::accu(arma::square(beta));
    }
    // Candidate neighborhoods already run in parallel. Avoid nested teams and
    // unused scratch arrays; a serialized nested region always has thread id0.
    const int workers=omp_in_parallel() ? 1 : omp_get_max_threads();
    std::vector<arma::vec> gradients(workers,arma::vec(P,arma::fill::zeros));
    std::vector<double> values(workers,0.0);
#pragma omp parallel for schedule(static) if(workers>1) num_threads(workers)
    for(int i=0;i<static_cast<int>(n);++i) {
      const int tid=omp_get_thread_num();
      const double own_eta=context->own_intercept+arma::dot(X_j.row(i),beta.col(0));
      const double lm=context->log_other(i),top=std::max(own_eta,lm);
      const double ls=top+std::log(std::exp(own_eta-top)+std::exp(lm-top));
      const DMRising own=dm_runtime_value_first(own_eta,context->own_count(i));
      const DMRising total=dm_runtime_value_first(ls,context->total(i));
      long double ll=-total.correction;
      if(context->own_count(i)>0)ll+=static_cast<long double>(context->own_count(i))*(own_eta-ls)+own.correction;
      ll+=context->fixed[i]+static_cast<long double>(context->other_count(i))*(lm-ls);
      const double row=static_cast<double>(ll);
      if(!std::isfinite(row)) {values[tid]=std::numeric_limits<double>::infinity();continue;}
      const double weight=std::exp(own_eta-ls);
      gradients[tid]-=(own.first-weight*total.first)*X_j.row(i).t();
      values[tid]-=row;
    }
    double value=0;
    for(int t=0;t<workers;++t) {gradient+=gradients[t];value+=values[t];}
    gradient+=lambda*beta;value+=.5*lambda*arma::accu(arma::square(beta));
    return value;
  }

  arma::mat ComputeHessian(const arma::mat& beta) {
    if(P==0)return arma::mat(0,0);
    if(context->single_category) {
      arma::vec eta=X_j*beta;eta+=context->own_intercept;
      if(!context->valid || !eta.is_finite()) {
        arma::mat bad(P,P);bad.fill(std::numeric_limits<double>::quiet_NaN());return bad;
      }
      return -lambda*arma::eye(P,P);
    }
    const int workers=omp_in_parallel() ? 1 : omp_get_max_threads();
    std::vector<arma::mat> blocks(workers,arma::mat(P,P,arma::fill::zeros));
#pragma omp parallel for schedule(static) if(workers>1) num_threads(workers)
    for(int i=0;i<static_cast<int>(n);++i) {
      const int tid=omp_get_thread_num();
      const double own_eta=context->own_intercept+arma::dot(X_j.row(i),beta.col(0));
      const double lm=context->log_other(i),top=std::max(own_eta,lm);
      const double ls=top+std::log(std::exp(own_eta-top)+std::exp(lm-top));
      const DMRising own=dm_runtime_derivatives(own_eta,context->own_count(i));
      const DMRising total=dm_runtime_derivatives(ls,context->total(i));
      const double weight=std::exp(own_eta-ls);
      const double curvature=own.second-weight*total.first+weight*weight*(total.first-total.second);
      blocks[tid]+=curvature*(X_j.row(i).t()*X_j.row(i));
    }
    arma::mat h(P,P,arma::fill::zeros);
    for(int t=0;t<workers;++t)h+=blocks[t];
    h.diag()-=lambda;return arma::symmatu(h);
  }
private:
  const arma::mat& X_j;
  const double lambda;
  const arma::uword n;
  const int P;
  const std::shared_ptr<const DMCachedConditionalContext> context;
};

// Facade keeps every downstream statistical operation identical. The reference
// implementation is retained verbatim for controlled side-by-side profiling.
class DMRegressionFunction {
public:
  DMRegressionFunction(const arma::mat& X,const arma::imat& Y,const arma::rowvec& alpha,
      const arma::ivec& PP,const arma::mat& Xbeta,int j,double lambda,const std::string* key=nullptr,
      std::shared_ptr<const DMCachedConditionalContext> shared_context=nullptr) {
    dm_runtime::Scope timer(dm_runtime::construction);
    if(j<0 || j>=static_cast<int>(Y.n_cols) || alpha.n_elem!=Y.n_cols || Xbeta.n_rows!=Y.n_rows ||
       Xbeta.n_cols!=Y.n_cols || X.n_rows!=Y.n_rows || PP.n_elem!=Y.n_cols ||
       PP(j)!=static_cast<int>(X.n_cols))Rcpp::stop("Invalid DM conditional dimensions in computational copy");
    if(dm_runtime::use_cached)cached.reset(new DMRegressionFunctionCached(X,Y,alpha,PP,Xbeta,j,lambda,key,shared_context));
    else reference.reset(new DMRegressionFunctionFullReference(X,Y,alpha,PP,Xbeta,j,lambda));
  }
  double EvaluateWithGradient(const arma::mat& beta,arma::mat& gradient) {
    dm_runtime::Scope timer(dm_runtime::objective);
    return cached ? cached->EvaluateWithGradient(beta,gradient) : reference->EvaluateWithGradient(beta,gradient);
  }
  arma::mat ComputeHessian(const arma::mat& beta) {
    dm_runtime::Scope timer(dm_runtime::hessian);
    return cached ? cached->ComputeHessian(beta) : reference->ComputeHessian(beta);
  }
private:
  std::unique_ptr<DMRegressionFunctionCached> cached;
  std::unique_ptr<DMRegressionFunctionFullReference> reference;
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
                       const int maxit = 10,
                       const bool state_start = false,
                       const std::string* supplied_context = nullptr,
                       std::shared_ptr<const DMCachedConditionalContext> shared_context=nullptr) {
  
  dm_runtime::Scope fit_timer(dm_runtime::fitting);
  
  // A truncated optimizer must define a reproducible proposal at each state.
  // Historical warm starts (and random cold starts) do not do that.  Memoize
  // completed deterministic fits only, keyed by their exact conditional inputs.
  // The category being fitted is excluded from the conditioning state.
  struct FitCache {
    std::string context;
    std::unordered_map<std::string, arma::vec> fits;
  };
  static thread_local FitCache cache;
  auto append = [](std::string& key, const void* data, size_t bytes) {
    if (bytes) key.append(static_cast<const char*>(data), bytes);
  };
  const std::string generated=supplied_context ? std::string() : dm_runtime_context_key(Y,intercepts,Xbeta,j);
  const std::string& context=supplied_context ? *supplied_context : generated;
  if (context != cache.context) {
    cache.context = context;
    cache.fits.clear();
  }
  std::string key;
  append(key, &lambda, sizeof(lambda));
  append(key, &maxit, sizeof(maxit));
  append(key, &state_start, sizeof(state_start));
  if (state_start) {
    if (beta_init.n_elem != X_j.n_cols || !beta_init.is_finite())
      Rcpp::stop("Invalid state-dependent coefficient start");
    append(key, beta_init.memptr(), beta_init.n_elem * sizeof(double));
  }
  const arma::uword design_shape[] = {X_j.n_rows, X_j.n_cols};
  append(key, design_shape, sizeof(design_shape));
  append(key, X_j.memptr(), X_j.n_elem * sizeof(double));
  auto cached = cache.fits.find(key);
  if (cached != cache.fits.end()) {
    if(dm_runtime::enabled)++dm_runtime::fit_cache_hits;
    return cached->second;
  }
  DMRegressionFunction dmrf(X_j,Y,intercepts,PP,Xbeta,j,lambda,&cache.context,shared_context);
  arma::vec beta_j = state_start ? beta_init : arma::vec(X_j.n_cols, arma::fill::zeros);
  
  ens::L_BFGS optimizer;
  optimizer.MaxIterations() = maxit;
  optimizer.Optimize(dmrf,beta_j);
  // L-BFGS may finish a failed line search with a nonfinite point.  Fall back
  // to its deterministic zero start, never to a history-dependent estimate.
  arma::mat gradient;
  if (!beta_j.is_finite() || !std::isfinite(dmrf.EvaluateWithGradient(beta_j, gradient)) ||
      !gradient.is_finite()) {
    ++dm_numeric_diagnostics::optimizer_fallbacks;
    beta_j = state_start ? beta_init : arma::vec(X_j.n_cols, arma::fill::zeros);
  }
  if (cache.fits.size() >= 256) cache.fits.clear();
  cache.fits.emplace(std::move(key), beta_j);
  
  return beta_j;
}

// Project coefficients from the SOURCE state onto the TARGET model. Shared
// predictors retain their values; births initialize at zero. Never use history.
arma::vec dm_project_state_start(const arma::uvec& source,
    const arma::vec& coefficients, const arma::uvec& target) {
  if (source.n_elem != coefficients.n_elem) Rcpp::stop("Invalid coefficient mapping");
  arma::vec start(target.n_elem, arma::fill::zeros);
  for (arma::uword k=0;k<target.n_elem;++k)
    for (arma::uword l=0;l<source.n_elem;++l)
      if (target(k)==source(l)) { start(k)=coefficients(l); break; }
  return start;
}

// [[Rcpp::export]]
arma::vec audit_dm_project_state_start(const arma::ivec& source,
    const arma::vec& coefficients, const arma::ivec& target) {
  if (arma::any(source<=0)||arma::any(target<=0)) Rcpp::stop("Indices must be one-based");
  return dm_project_state_start(arma::conv_to<arma::uvec>::from(source-1),
    coefficients,arma::conv_to<arma::uvec>::from(target-1));
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

std::string vectorToKey(const arma::ivec& config);

struct DMLaplaceFit {
  arma::vec mode;
  arma::mat precision;
  double log_score;
};

double dm_indicator_log_prior(const arma::ivec& xi_j,
                              double a, double b) {
  return arma::accu(
    arma::lgamma(arma::conv_to<arma::vec>::from(xi_j) + a) +
      arma::lgamma(1.0 - arma::conv_to<arma::vec>::from(xi_j) + b)
  );
}

// [[Rcpp::export]]
double audit_dm_indicator_log_prior(const arma::ivec& xi,double a=1,double b=9,
    bool legacy_integer_arithmetic=false) {
  // Legacy branch is diagnostic-only, to reproduce the original truncation.
  if (legacy_integer_arithmetic)
    return arma::sum(arma::lgamma(xi+a)+arma::lgamma(1.0-xi+b));
  return dm_indicator_log_prior(xi,a,b);
}

arma::mat stabilize_dm_precision(const arma::mat& precision_raw,
                                 double& log_determinant) {
  if (precision_raw.n_rows == 0) {
    log_determinant = 0.0;
    return precision_raw;
  }
  arma::mat precision = arma::symmatu(precision_raw);
  if (!precision.is_finite()) Rcpp::stop("Nonfinite Laplace posterior precision");
  arma::mat chol_factor;
  if (arma::chol(chol_factor, precision)) {
    log_determinant = 2.0 * arma::accu(arma::log(chol_factor.diag()));
    return precision;
  }
  double jitter = 1e-8;
  for (int attempt = 0; attempt < 9; ++attempt) {
    precision = arma::symmatu(precision_raw);
    precision.diag() += jitter;
    if (arma::chol(chol_factor, precision)) {
      log_determinant = 2.0 * arma::accu(arma::log(chol_factor.diag()));
      return precision;
    }
    jitter *= 10.0;
  }
  Rcpp::stop("Laplace posterior precision is not positive definite after stabilization.");
  return precision;
}

DMLaplaceFit fit_dm_conditional_laplace(
    const arma::ivec& xi_j, int j, const arma::ivec& PP,
    const arma::mat& Xbeta, const arma::imat& Y, const arma::mat& X,
    const arma::vec& alpha, const arma::mat& alpha_mat,
    const arma::ivec& y_sum_i, double log_lik_current,
    const arma::mat& gamma_current, double a, double b, double r2_j,
    int optimizer_maxit, const arma::vec& warm_start) {
  const int n = Y.n_rows;
  arma::uvec selected = arma::find(xi_j > 0);
  arma::ivec PP_new = PP;
  PP_new(j) = static_cast<int>(selected.n_elem);
  arma::mat X_new = X.cols(selected);
  arma::vec mode;
  arma::mat precision;
  double log_determinant = 0.0;
  const double posterior_precision = 1.0 / r2_j;

  if (selected.n_elem > 0) {
    arma::vec start = warm_start.n_elem == selected.n_elem ?
      warm_start : arma::vec(selected.n_elem, arma::fill::zeros);
    mode = dm_reg_lbfgs(
      X_new, Y, alpha.t(), PP_new, Xbeta, j, posterior_precision,
      start, optimizer_maxit
    );
    arma::mat hessian = dm_reg_hessian(
      X_new, Y, alpha.t(), PP_new, Xbeta, j, mode, posterior_precision
    );
    precision = stabilize_dm_precision(-hessian, log_determinant);
  } else {
    mode = arma::vec();
    precision = arma::mat(0, 0);
  }

  arma::vec replacement(n,arma::fill::zeros);
  if (selected.n_elem>0) replacement=X_new*mode;
  double log_likelihood_at_mode=log_lik_current+
    dm_stable_replacement_delta(Y,alpha_mat,Xbeta,j,replacement);
  const double dimension = static_cast<double>(selected.n_elem);
  double log_beta_prior = -0.5 * dimension *
    (std::log(2.0 * M_PI) + std::log(r2_j));
  if (selected.n_elem > 0) {
    log_beta_prior -= 0.5 * arma::dot(mode, mode) / r2_j;
  }
  double log_laplace_correction =
    0.5 * dimension * std::log(2.0 * M_PI) - 0.5 * log_determinant;
  double log_score = dm_indicator_log_prior(xi_j, a, b) +
    log_likelihood_at_mode + log_beta_prior + log_laplace_correction;
  return DMLaplaceFit{mode, precision, log_score};
}

double dm_log_barker(double log_ratio) {
  if (log_ratio > 0.0) return -std::log1p(std::exp(-log_ratio));
  return log_ratio - std::log1p(std::exp(log_ratio));
}

arma::vec proposal_gagnon_laplace(
    const arma::ivec& xi_j, int j, const arma::ivec& PP,
    const arma::mat& Xbeta,
    std::unordered_map<std::string, arma::vec>& beta_j_hat_memory,
    const arma::imat& Y, const arma::mat& X, const arma::vec& alpha,
    const arma::mat& alpha_mat, const arma::ivec& y_sum_i,
    double log_lik_current, const arma::mat& gamma_current,
    double a, double b, double r2_j, int optimizer_maxit,
    bool cache_fits) {
  const int P = xi_j.n_elem;
  std::string current_key = vectorToKey(xi_j);
  arma::vec current_warm;
  auto current_it = beta_j_hat_memory.find(current_key);
  if (current_it != beta_j_hat_memory.end()) current_warm = current_it->second;
  DMLaplaceFit current_fit = fit_dm_conditional_laplace(
    xi_j, j, PP, Xbeta, Y, X, alpha, alpha_mat, y_sum_i,
    log_lik_current, gamma_current, a, b, r2_j, optimizer_maxit,
    current_warm
  );

  std::vector<std::string> candidate_keys(P);
  std::vector<arma::vec> candidate_warm(P);
  for (int p = 0; p < P; ++p) {
    arma::ivec candidate = xi_j;
    candidate(p) = 1 - candidate(p);
    candidate_keys[p] = vectorToKey(candidate);
    auto it = beta_j_hat_memory.find(candidate_keys[p]);
    if (it != beta_j_hat_memory.end()) candidate_warm[p] = it->second;
  }

  arma::vec log_weights(P, arma::fill::zeros);
  std::vector<arma::vec> candidate_modes(P);
#pragma omp parallel for schedule(dynamic)
  for (int p = 0; p < P; ++p) {
    arma::ivec candidate = xi_j;
    candidate(p) = 1 - candidate(p);
    DMLaplaceFit candidate_fit = fit_dm_conditional_laplace(
      candidate, j, PP, Xbeta, Y, X, alpha, alpha_mat, y_sum_i,
      log_lik_current, gamma_current, a, b, r2_j, optimizer_maxit,
      candidate_warm[p]
    );
    candidate_modes[p] = candidate_fit.mode;
    log_weights(p) = dm_log_barker(candidate_fit.log_score - current_fit.log_score);
  }

  if (cache_fits) {
    beta_j_hat_memory[current_key] = current_fit.mode;
    for (int p = 0; p < P; ++p) {
      beta_j_hat_memory[candidate_keys[p]] = candidate_modes[p];
    }
  }
  double maximum = log_weights.max();
  arma::vec probabilities = arma::exp(log_weights - maximum);
  probabilities /= arma::accu(probabilities);
  return probabilities;
}

// Deterministic statewise audit entry point for the Gagnon-adapted model
// proposal.  beta_full is P by J with zeros in inactive positions.
// [[Rcpp::export]]
arma::mat audit_dm_gagnon_laplace_probabilities(
    const arma::mat& X, const arma::imat& Y, const arma::vec& alpha,
    const arma::imat& xi, const arma::mat& beta_full,
    const arma::vec& r2, double a = 1.0, double b = 9.0,
    int optimizer_maxit = 25) {
  const int n = Y.n_rows;
  const int P = X.n_cols;
  const int J = Y.n_cols;
  if (xi.n_rows != static_cast<arma::uword>(P) ||
      xi.n_cols != static_cast<arma::uword>(J) ||
      beta_full.n_rows != static_cast<arma::uword>(P) ||
      beta_full.n_cols != static_cast<arma::uword>(J)) {
    Rcpp::stop("xi and beta_full must both have dimensions P by J.");
  }
  if (alpha.n_elem != static_cast<arma::uword>(J) ||
      r2.n_elem != static_cast<arma::uword>(J)) {
    Rcpp::stop("alpha and r2 must have length J.");
  }
  if (optimizer_maxit <= 0) Rcpp::stop("optimizer_maxit must be positive.");
  arma::mat masked_beta = beta_full % arma::conv_to<arma::mat>::from(xi);
  arma::mat Xbeta = X * masked_beta;
  arma::ivec PP = arma::sum(xi, 0).t();
  arma::mat alpha_mat(n, J);
  alpha_mat.each_row() = alpha.t();
  arma::mat gamma = arma::trunc_exp(alpha_mat + Xbeta);
  arma::ivec y_sum_i = arma::sum(Y, 1);
  arma::vec gamma_sum_i = arma::sum(gamma, 1);
  double log_lik=dm_stable_log_kernel(Y,alpha_mat+Xbeta);
  arma::mat probabilities(P, J, arma::fill::zeros);
  for (int j = 0; j < J; ++j) {
    std::unordered_map<std::string, arma::vec> memory;
    arma::uvec selected = arma::find(xi.col(j) > 0);
    if (selected.n_elem > 0) {
      arma::vec beta_column = masked_beta.col(j);
      memory[vectorToKey(xi.col(j))] = beta_column.elem(selected);
    }
    probabilities.col(j) = proposal_gagnon_laplace(
      xi.col(j), j, PP, Xbeta, memory, Y, X, alpha, alpha_mat,
      y_sum_i, log_lik, gamma, a, b, r2(j), optimizer_maxit, false
    );
  }
  return probabilities;
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

double dmvn_prec(const arma::vec& x,const arma::vec& mean,
                  const arma::mat& Q,bool log_result=true) {
  if (!x.is_finite()||!mean.is_finite())
    return log_result ? -std::numeric_limits<double>::infinity() : 0.0;
  arma::mat precision=dm_gaussian_precision(Q),L;
  arma::chol(L,precision,"lower");
  arma::vec z=L.t()*(x-mean);
  double value=-0.5*(x.n_elem*std::log(2.0*M_PI)
    -2.0*arma::accu(arma::log(L.diag()))+arma::dot(z,z));
  return log_result ? value : std::exp(value);
}

arma::vec rmvn_prec(const arma::vec& mu,const arma::mat& Q) {
  if (!mu.is_finite()) Rcpp::stop("Nonfinite DM Gaussian proposal mean");
  arma::mat precision=dm_gaussian_precision(Q),L;
  arma::chol(L,precision,"lower");
  return mu+arma::solve(arma::trimatu(L.t()),arma::randn(mu.n_elem));
}

struct DMWithinTerms {
  double likelihood_delta, prior_delta, log_forward, log_reverse, ratio;
};

DMWithinTerms dm_within_terms(const arma::mat& design,const arma::imat& Y,
    const arma::vec& alpha,const arma::ivec& pp,const arma::mat& xb,int j,
    const arma::vec& current,const arma::vec& candidate,double r2,double scale,
    const arma::mat& qf, bool fixed_design=false) {
  arma::mat next_xb=xb; next_xb.col(j)=design*candidate;
  arma::mat qr=qf;
  if (!fixed_design) qr=dm_gaussian_precision(-dm_reg_hessian(
    design,Y,alpha.t(),pp,next_xb,j,candidate,1.0/r2))/(scale*scale);
  arma::mat alpha_mat(Y.n_rows,Y.n_cols); alpha_mat.each_row()=alpha.t();
  DMWithinTerms terms;
  terms.likelihood_delta=dm_stable_replacement_delta(Y,alpha_mat,xb,j,next_xb.col(j));
  terms.prior_delta=-0.5*(arma::dot(candidate,candidate)-arma::dot(current,current))/r2;
  terms.log_forward=dmvn_prec(candidate,current,qf);
  terms.log_reverse=dmvn_prec(current,candidate,qr);
  terms.ratio=terms.likelihood_delta+terms.prior_delta+
    (fixed_design ? 0.0 : terms.log_reverse-terms.log_forward);
  return terms;
}

// Fixed-pair audit of the SAME MH calculation used by production; no RNG.
// [[Rcpp::export]]
Rcpp::List audit_dm_within_model_terms(const arma::mat& design,const arma::imat& Y,
    const arma::vec& alpha,const arma::mat& xb,int category,const arma::vec& current,
    const arma::vec& candidate,double r2=10,double scale=0.25,bool fixed_design=false) {
  int j=category-1;
  if (j<0||j>=static_cast<int>(Y.n_cols)||design.n_cols!=current.n_elem||
      current.n_elem!=candidate.n_elem||current.n_elem==0||r2<=0||scale<=0)
    Rcpp::stop("Invalid within-model audit inputs");
  arma::ivec pp(Y.n_cols,arma::fill::zeros); pp(j)=current.n_elem;
  arma::mat raw= fixed_design ? arma::mat(design.t()*design/design.n_rows+
    arma::eye(design.n_cols,design.n_cols)/r2) : arma::mat(-dm_reg_hessian(
    design,Y,alpha.t(),pp,xb,j,current,1.0/r2));
  arma::mat qf=dm_gaussian_precision(raw)/(scale*scale);
  DMWithinTerms terms=dm_within_terms(design,Y,alpha,pp,xb,j,current,candidate,r2,scale,qf,fixed_design);
  return Rcpp::List::create(Rcpp::Named("likelihood_delta")=terms.likelihood_delta,
    Rcpp::Named("prior_delta")=terms.prior_delta,Rcpp::Named("log_forward")=terms.log_forward,
    Rcpp::Named("log_reverse")=terms.log_reverse,Rcpp::Named("ratio")=terms.ratio,
    Rcpp::Named("forward_precision")=qf);
}

#include "joint_kernel.h"

// Fixed-pair entry point: no optimizer, RNG or chain advancement.
// [[Rcpp::export]]
Rcpp::List audit_dm_joint_block_terms(const arma::mat& design,const arma::imat& Y,
    const arma::vec& alpha,const arma::mat& xb,int category,const arma::vec& beta,
    const arma::vec& candidate,double s2=10,double r2=10,double scale=1.0) {
  return dm_joint::audit(design,Y,alpha,xb,category,beta,candidate,s2,r2,scale);
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

arma::vec normalize_dm_similarity_weights(const arma::vec& log_weights,
                                          const arma::ivec& indicators,
                                          bool balance_directions,
                                          double uniform_mixture) {
  int P = indicators.n_elem;
  if (!log_weights.is_finite()) Rcpp::stop("Nonfinite DM similarity weights");
  arma::vec probabilities(P, arma::fill::zeros);
  if (!balance_directions) {
    probabilities = arma::normalise(
      arma::trunc_exp(log_weights - log_weights.max()), 1
    );
  } else {
    arma::uvec births = arma::find(indicators == 0);
    arma::uvec deaths = arma::find(indicators == 1);
    bool has_births = births.n_elem > 0;
    bool has_deaths = deaths.n_elem > 0;
    double birth_mass = has_births && has_deaths ? 0.5 : (has_births ? 1.0 : 0.0);
    double death_mass = has_births && has_deaths ? 0.5 : (has_deaths ? 1.0 : 0.0);
    if (has_births) {
      arma::vec z = arma::trunc_exp(
        log_weights.elem(births) - log_weights.elem(births).max()
      );
      probabilities.elem(births) = birth_mass * z / arma::accu(z);
    }
    if (has_deaths) {
      arma::vec z = arma::trunc_exp(
        log_weights.elem(deaths) - log_weights.elem(deaths).max()
      );
      probabilities.elem(deaths) = death_mass * z / arma::accu(z);
    }
  }

  if (uniform_mixture > 0.0) {
    probabilities *= (1.0 - uniform_mixture);
    if (balance_directions) {
      arma::uvec births = arma::find(indicators == 0);
      arma::uvec deaths = arma::find(indicators == 1);
      bool has_births = births.n_elem > 0;
      bool has_deaths = deaths.n_elem > 0;
      double birth_mass = has_births && has_deaths ? 0.5 : (has_births ? 1.0 : 0.0);
      double death_mass = has_births && has_deaths ? 0.5 : (has_deaths ? 1.0 : 0.0);
      if (has_births) {
        probabilities.elem(births) += birth_mass * uniform_mixture /
          static_cast<double>(births.n_elem);
      }
      if (has_deaths) {
        probabilities.elem(deaths) += death_mass * uniform_mixture /
          static_cast<double>(deaths.n_elem);
      }
    } else {
      probabilities += uniform_mixture / static_cast<double>(P);
    }
  }
  return probabilities;
}

// [[Rcpp::export]]
arma::vec audit_dm_similarity_normalization(const arma::vec& log_weights,
                                            const arma::ivec& indicators,
                                            bool balance_directions = true,
                                            double uniform_mixture = 0.02) {
  if (log_weights.n_elem != indicators.n_elem || log_weights.n_elem == 0) {
    Rcpp::stop("log_weights and indicators must have the same positive length.");
  }
  if (arma::any(indicators < 0) || arma::any(indicators > 1)) {
    Rcpp::stop("indicators must contain only zeros and ones.");
  }
  if (uniform_mixture < 0.0 || uniform_mixture >= 1.0) {
    Rcpp::stop("uniform_mixture must lie in [0, 1).");
  }
  return normalize_dm_similarity_weights(
    log_weights, indicators, balance_directions, uniform_mixture
  );
}

arma::vec proposal_inf(const arma::ivec& xi_j, const int j,
                       const arma::ivec& PP, const arma::mat& Xbeta,
                       std::unordered_map<std::string, arma::vec>& beta_j_hat_memory,
                       const arma::imat& Y, const arma::mat& X,
                       const arma::vec& alpha, const arma::mat& alpha_mat,
                       const arma::ivec& y_sum_i, 
                       const double log_lik_current,  
                       const arma::mat& gamma_current, 
                       double g, double c,
                       bool balance_directions,
                       double uniform_mixture,
                       int score_optimizer_maxit,
                       bool cache_score_fits) {
  dm_runtime::Scope neighborhood_timer(dm_runtime::neighborhood);
  
  int P = xi_j.n_elem;
  int n = Y.n_rows;
  // double h = 0.01;
  
  arma::ivec Yj = Y.col(j);
  
  const double log_lik0=dm_stable_log_kernel(Y,alpha_mat);
  std::unique_ptr<const DMRuntimeReplacementContext> replacement_context;
  if(dm_runtime::use_cached && dm_runtime::cache_replacements)
    replacement_context.reset(new DMRuntimeReplacementContext(Y,alpha_mat,Xbeta,j));
  const bool share_context=dm_runtime::use_cached && dm_runtime::share_neighborhood_context;
  const std::string context_key=share_context ? dm_runtime_context_key(Y,alpha.t(),Xbeta,j) : std::string();
  std::shared_ptr<const DMCachedConditionalContext> conditional_context;
  if(share_context) {
    conditional_context=std::make_shared<DMCachedConditionalContext>(Y,alpha.t(),Xbeta,j);
    if(dm_runtime::enabled)++dm_runtime::conditional_context_misses;
  }
  // ---- Initialize proposal vectors ----
  arma::vec log_prop(P, arma::fill::zeros);
  std::vector<std::string> candidate_configs(P);
  std::vector<arma::vec> candidate_beta_hats(P);
  std::vector<int> candidate_has_beta(P, 0);
  std::vector<int> invalid_candidate(P, 0);
  
  // ---- Parallel loop ----
#pragma omp parallel for schedule(dynamic)
  for (int p = 0; p < P; p++) {
    
    arma::ivec xi_j_new = xi_j;
    xi_j_new(p) = 1 - xi_j(p);
    arma::uvec selected_new = arma::find(xi_j_new > 0);
    std::string config_new = vectorToKey(arma::conv_to<arma::ivec>::from(xi_j_new));
    candidate_configs[p] = config_new;
    
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
        beta_hat_j = dm_reg_lbfgs(X_new, Y, alpha.t(), PP_new, Xbeta, j, lambda_new, it->second, score_optimizer_maxit,
          false,share_context ? &context_key : nullptr,conditional_context);
      } else {
        // New configuration: start from random
        arma::vec beta_start(PP_new(j), arma::fill::zeros);
        beta_hat_j = dm_reg_lbfgs(X_new, Y, alpha.t(), PP_new, Xbeta, j, lambda_new, beta_start, score_optimizer_maxit,
          false,share_context ? &context_key : nullptr,conditional_context);
      }
      candidate_beta_hats[p] = beta_hat_j;
      candidate_has_beta[p] = 1;
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
    double log_lik_proposed=log_lik_current+
      (replacement_context ? replacement_context->delta(Xbeta_j_new) :
        dm_stable_replacement_delta(Y,alpha_mat,Xbeta,j,Xbeta_j_new));

    // ---- Likelihood ratio statistic (proposed vs intercept-only) ----
    double lik_ratio_stat = -2.0 * (log_lik0 - log_lik_proposed);
    if (!std::isfinite(lik_ratio_stat)) {
      invalid_candidate[p] = 1;
      continue;
    }
    if (lik_ratio_stat < 0.0)
      lik_ratio_stat = 0.0;
    
    // ---- Compute -log10(p-value) ----
    double log_pval = selected_new.n_elem > 0 ?
      R::pchisq(lik_ratio_stat, selected_new.n_elem, false, true) : 0.0;
    double log10_LR_pvalue = -log_pval / std::log(10.0);
    
    // ---- Proposal weight ----
    log_prop(p) = arma::trunc_exp(g * arma::trunc_log(log10_LR_pvalue));
  }

  if (std::any_of(invalid_candidate.begin(),invalid_candidate.end(),
                  [](int value) { return value != 0; }))
    Rcpp::stop("Nonfinite DM flip candidate likelihood");
  // Cache writes are deliberately performed after the parallel region.  The
  // candidate fits are independent, but concurrent unordered_map mutation is
  // not safe. Historical entries are retained for API compatibility, but the
  // optimizer uses only its exact-context cache, never these as warm starts.
  if (cache_score_fits) {
    for (int p = 0; p < P; ++p) {
      if (candidate_has_beta[p]) {
        beta_j_hat_memory[candidate_configs[p]] = candidate_beta_hats[p];
      }
    }
  }
  
  return normalize_dm_similarity_weights(
    log_prop, xi_j, balance_directions, uniform_mixture
  );
}

double intercept_only_loglik(const arma::imat& Y,
                             const arma::mat& alpha_mat,
                             const arma::ivec& y_sum_i) {
  return dm_stable_log_kernel(Y,alpha_mat);
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
    double c,
    int score_optimizer_maxit) {
  dm_runtime::Scope neighborhood_timer(dm_runtime::neighborhood);
  int m = static_cast<int>(inactive_neighbors.size());
  if (m == 0) {
    return arma::vec();
  }
  
  int n = Y.n_rows;
  arma::ivec Yj = Y.col(j);
  double log_lik0 = intercept_only_loglik(Y, alpha_mat, y_sum_i);
  std::unique_ptr<const DMRuntimeReplacementContext> replacement_context;
  if(dm_runtime::use_cached && dm_runtime::cache_replacements)
    replacement_context.reset(new DMRuntimeReplacementContext(Y,alpha_mat,Xbeta,j));
  const bool share_context=dm_runtime::use_cached && dm_runtime::share_neighborhood_context;
  const std::string context_key=share_context ? dm_runtime_context_key(Y,alpha.t(),Xbeta,j) : std::string();
  std::shared_ptr<const DMCachedConditionalContext> conditional_context;
  if(share_context) {
    conditional_context=std::make_shared<DMCachedConditionalContext>(Y,alpha.t(),Xbeta,j);
    if(dm_runtime::enabled)++dm_runtime::conditional_context_misses;
  }
  
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
          X_move, Y, alpha.t(), PP_move, Xbeta, j, lambda_move_j, it->second, score_optimizer_maxit,
          false,share_context ? &context_key : nullptr,conditional_context
        );
      } else {
        arma::vec beta_start(PP_move(j), arma::fill::zeros);
        beta_hat_move = dm_reg_lbfgs(
          X_move, Y, alpha.t(), PP_move, Xbeta, j, lambda_move_j, beta_start, score_optimizer_maxit,
          false,share_context ? &context_key : nullptr,conditional_context
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
    
    double log_lik_move=log_lik_current+
      (replacement_context ? replacement_context->delta(Xbeta_j_move) :
        dm_stable_replacement_delta(Y,alpha_mat,Xbeta,j,Xbeta_j_move));
    double lik_ratio_stat = -2.0 * (log_lik0 - log_lik_move);
    if (!std::isfinite(lik_ratio_stat))
      Rcpp::stop("Nonfinite DM local candidate likelihood");
    if (lik_ratio_stat < 0.0) {
      lik_ratio_stat = 0.0;
    }
    
    int df = std::max(1, static_cast<int>(selected_move.n_elem));
    double log_pval = R::pchisq(lik_ratio_stat, df, false, true);
    double log10_lrt_pvalue = -log_pval / std::log(10.0);
    score(idx) = arma::trunc_exp(lambda_move * arma::trunc_log(log10_lrt_pvalue));
  }
  
  if (!score.is_finite()) Rcpp::stop("Nonfinite DM local similarity weights");
  return arma::normalise(arma::trunc_exp(score - score.max()), 1);
}

// Numerical audit entry points: these call the production calculations.
// [[Rcpp::export]]
arma::mat audit_dm_runtime_v2_rising(const arma::vec& log_a,const arma::ivec& counts) {
  if(log_a.n_elem!=counts.n_elem)Rcpp::stop("audit lengths differ");
  arma::mat answer(counts.n_elem,6);
  for(arma::uword i=0;i<counts.n_elem;++i) {
    DMRising vf=dm_rising_value_first(log_a(i),counts(i));
    DMRising d=dm_rising_derivatives(log_a(i),counts(i));
    DMRising full=dm_rising(log_a(i),counts(i));
    answer(i,0)=dm_rising_value(log_a(i),counts(i));
    answer(i,1)=vf.correction;answer(i,2)=vf.first;
    answer(i,3)=d.first;answer(i,4)=d.second;answer(i,5)=full.second;
  }
  return answer;
}

// [[Rcpp::export]]
Rcpp::List audit_dm_runtime_v2_replacement(const arma::imat& Y,const arma::mat& alpha_mat,
    const arma::mat& xb,int category,const arma::mat& replacements) {
  const int j=category-1;
  if(j<0||j>=static_cast<int>(Y.n_cols)||Y.n_rows!=alpha_mat.n_rows||Y.n_cols!=alpha_mat.n_cols||
      Y.n_rows!=xb.n_rows||Y.n_cols!=xb.n_cols||Y.n_rows!=replacements.n_rows||Y.min()<0)
    Rcpp::stop("invalid replacement audit dimensions/counts");
  DMRuntimeReplacementContext context(Y,alpha_mat,xb,j);
  arma::vec cached(replacements.n_cols),reference(replacements.n_cols);
  for(arma::uword k=0;k<replacements.n_cols;++k) {
    cached(k)=context.delta(replacements.col(k));
    reference(k)=dm_stable_replacement_delta(Y,alpha_mat,xb,j,replacements.col(k));
  }
  return Rcpp::List::create(Rcpp::Named("cached")=cached,Rcpp::Named("reference")=reference);
}

// [[Rcpp::export]]
arma::mat audit_dm_rising_values(const arma::vec& log_a, const arma::ivec& counts) {
  if (log_a.n_elem != counts.n_elem) Rcpp::stop("audit lengths differ");
  arma::mat values(counts.n_elem,3);
  for (arma::uword i=0;i<counts.n_elem;++i) {
    DMRising value=dm_rising(log_a(i),counts(i));
    values(i,0)=value.correction;
    values(i,1)=value.first;
    values(i,2)=value.second;
  }
  return values;
}

// [[Rcpp::export]]
double audit_dm_stable_likelihood(const arma::imat& Y,const arma::mat& eta) {
  if (Y.n_rows!=eta.n_rows||Y.n_cols!=eta.n_cols||Y.min()<0)
    Rcpp::stop("invalid audit data");
  return dm_stable_log_kernel(Y,eta);
}

// [[Rcpp::export]]
Rcpp::List audit_dm_conditional_numerics(const arma::mat& X,
    const arma::imat& Y,const arma::vec& alpha,const arma::mat& Xbeta,
    int category,const arma::vec& beta,double lambda,int maxit=3,
    bool run_optimizer=false, bool state_start=false) {
  int j=category-1;
  if (j<0||j>=static_cast<int>(Y.n_cols)||beta.n_elem!=X.n_cols)
    Rcpp::stop("invalid conditional audit arguments");
  arma::ivec pp(Y.n_cols,arma::fill::zeros); pp(j)=X.n_cols;
  DMRegressionFunction objective(X,Y,alpha.t(),pp,Xbeta,j,lambda);
  arma::mat gradient;
  double value=objective.EvaluateWithGradient(beta,gradient);
  arma::vec mode=beta;
  if (run_optimizer) mode=dm_reg_lbfgs(X,Y,alpha.t(),pp,Xbeta,j,lambda,beta,maxit,state_start);
  arma::mat alpha_mat(Y.n_rows,Y.n_cols); alpha_mat.each_row()=alpha.t();
  return Rcpp::List::create(
    Rcpp::Named("value")=value,Rcpp::Named("gradient")=gradient,
    Rcpp::Named("hessian")=objective.ComputeHessian(beta),
    Rcpp::Named("mode")=mode,
    Rcpp::Named("mode_hessian")=objective.ComputeHessian(mode),
    Rcpp::Named("delta")=dm_stable_replacement_delta(Y,alpha_mat,Xbeta,j,X*beta));
}

// [[Rcpp::export]]
Rcpp::List audit_dm_gaussian_numerics(const arma::vec& mean,
    const arma::mat& precision,const arma::vec& x,int draws=0) {
  arma::mat samples(mean.n_elem,draws);
  for (int i=0;i<draws;++i) samples.col(i)=rmvn_prec(mean,precision);
  return Rcpp::List::create(
    Rcpp::Named("precision")=dm_gaussian_precision(precision),
    Rcpp::Named("log_density")=dmvn_prec(x,mean,precision),
    Rcpp::Named("samples")=samples);
}

// [[Rcpp::export]]
arma::vec audit_dm_safe_acceptance(const arma::vec& ratios) {
  arma::vec result(ratios.n_elem);
  for (arma::uword i=0;i<ratios.n_elem;++i) result(i)=dm_safe_log_acceptance(ratios(i));
  return result;
}

// [[Rcpp::export]]
Rcpp::List audit_dm_neighborhood_probabilities(const arma::mat& X,
    const arma::imat& Y,const arma::vec& alpha,const arma::imat& xi,
    const arma::mat& beta_full,int category,const arma::imat& adjacency,
    double g=1,double c=1,double mixture=0.02,int maxit=1) {
  int j=category-1, P=X.n_cols;
  arma::mat Xbeta=X*(beta_full%arma::conv_to<arma::mat>::from(xi));
  arma::mat alpha_mat(Y.n_rows,Y.n_cols); alpha_mat.each_row()=alpha.t();
  arma::ivec pp=arma::sum(xi,0).t(), totals=arma::sum(Y,1);
  arma::mat unused_gamma;
  double ll=dm_stable_log_kernel(Y,alpha_mat+Xbeta);
  std::unordered_map<std::string,arma::vec> memory;
  arma::vec flip=proposal_inf(xi.col(j),j,pp,Xbeta,memory,Y,X,alpha,
    alpha_mat,totals,ll,unused_gamma,g,c,true,mixture,maxit,true);
  std::vector<std::vector<int>> neighbors(P);
  for (int p=0;p<P;++p) for (int q=0;q<P;++q)
    if (p!=q&&(adjacency(p,q)||adjacency(q,p))) neighbors[p].push_back(q);
  std::vector<int> active=active_with_inactive_neighbors(xi.col(j),neighbors);
  arma::mat swaps(P,P,arma::fill::zeros);
  for (int drop:active) {
    std::vector<int> inactive=inactive_neighbors_for_active(drop,xi.col(j),neighbors);
    arma::vec weights=proposal_move_weights(xi.col(j),drop,inactive,j,pp,
      Xbeta,memory,Y,X,alpha,alpha_mat,totals,ll,unused_gamma,g,c,maxit);
    for (size_t k=0;k<inactive.size();++k)
      swaps(drop,inactive[k])=weights(k)/active.size();
  }
  return Rcpp::List::create(Rcpp::Named("flip")=flip,Rcpp::Named("swap")=swaps);
}

// [[Rcpp::export]]
Rcpp::List DM_joint_block_VS_RJMCMC(const arma::mat& X, const arma::imat& Y,
                                 const arma::vec s2, const arma::vec r2,
                                 double a, double b, double g_initial, double c,
                                 int runs, int burn, int thin = 1, 
                                 bool informed = true,
                                 double trgt_acpt_alpha = 0.234,
                                 bool adapt_g = true,
                                 int adaptation_start = 100,
                                 int adaptation_end = -1,
                                 int window_size = 25,
                                 double rm_scale = 0.5,
                                 double decay_power = 0.6,
                                 double min_g = 0.05,
                                 double max_g = 10.0,
                                 Rcpp::Nullable<Rcpp::List> neigh = R_NilValue,
                                 bool use_move = false,
                                 int n_move_attempts = 1,
                                 double lambda_move = 1.0,
                                 std::string adaptation_metric = "diversity",
                                 double target_effective_neighbors = -1.0,
                                 double adaptation_offset = 100.0,
                                 bool balance_directions = true,
                                 double uniform_mixture = 0.02,
                                 int score_optimizer_maxit = 1,
                                 int coefficient_optimizer_maxit = 3,
                                 Rcpp::Nullable<Rcpp::IntegerMatrix> initial_xi = R_NilValue,
                                 bool cache_score_fits = true,
                                 std::string discrete_proposal = "similarity",
                                 int laplace_optimizer_maxit = 25,
                                 bool coefficient_warm_start = false,
                                 bool within_model_update = false,
                                 int within_model_attempts = 1,
                                 double within_model_scale = 0.25,
                                 double within_model_fixed_probability = 0.0,
                                 double within_model_fixed_scale = 0.25,
                                 int stop_after = -1,
                                 Rcpp::Nullable<Rcpp::List> resume_state = R_NilValue,
                                 Rcpp::Nullable<Rcpp::List> resume_draws = R_NilValue,
                                 bool joint_update = false,
                                 double joint_scale = 1.0,
                                 int joint_attempts = 1,
                                 Rcpp::Nullable<Rcpp::List> pilot_start = R_NilValue){

  if(resume_state.isNotNull()||resume_draws.isNotNull()||stop_after>=0)
    Rcpp::stop("Experimental bounded pilot: use pilot_start, not production checkpoint continuation");
  if(joint_attempts<1||!std::isfinite(joint_scale)||joint_scale<=0)
    Rcpp::stop("Joint attempts/scale must be positive");
  if(pilot_start.isNotNull()&&(burn!=0||adapt_g))
    Rcpp::stop("Saved-endpoint pilots require burn=0 and frozen g/intercept adaptation");
  
  if (runs <= 0) {
    Rcpp::stop("runs must be positive.");
  }
  if (within_model_attempts < 1 || !std::isfinite(within_model_scale) || within_model_scale <= 0 ||
      !std::isfinite(within_model_fixed_probability) || within_model_fixed_probability<0 ||
      within_model_fixed_probability>1 || !std::isfinite(within_model_fixed_scale) || within_model_fixed_scale<=0)
    Rcpp::stop("within-model attempts and scale must be positive");
  if (coefficient_warm_start && discrete_proposal != "similarity")
    Rcpp::stop("State-dependent coefficient starts currently support similarity proposals only");
  if (burn < 0 || burn > runs) {
    Rcpp::stop("burn must lie between 0 and runs.");
  }
  if (thin < 0) {
    Rcpp::stop("thin must be non-negative.");
  }
  if (adaptation_metric != "diversity" && adaptation_metric != "legacy") {
    Rcpp::stop("adaptation_metric must be 'diversity' or 'legacy'.");
  }
  if (uniform_mixture < 0.0 || uniform_mixture >= 1.0) {
    Rcpp::stop("uniform_mixture must lie in [0, 1).");
  }
  if (adaptation_offset < 0.0) {
    Rcpp::stop("adaptation_offset must be non-negative.");
  }
  if (score_optimizer_maxit <= 0 || coefficient_optimizer_maxit <= 0) {
    Rcpp::stop("optimizer iteration limits must be positive.");
  }
  if (discrete_proposal != "similarity" &&
      discrete_proposal != "gagnon_laplace") {
    Rcpp::stop("discrete_proposal must be 'similarity' or 'gagnon_laplace'.");
  }
  if (laplace_optimizer_maxit <= 0) {
    Rcpp::stop("laplace_optimizer_maxit must be positive.");
  }
  if (discrete_proposal == "gagnon_laplace") {
    if (!informed) {
      Rcpp::stop("gagnon_laplace requires informed = true.");
    }
    if (adapt_g) {
      Rcpp::stop("gagnon_laplace has no similarity exponent; set adapt_g = false.");
    }
    if (use_move) {
      Rcpp::stop("gagnon_laplace currently implements the one-flip neighborhood only.");
    }
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
  // Diversity adaptation is deliberately online: by default it continues
  // through the retained run with diminishing gains.  Legacy hill climbing
  // retains the old burn-in-only default.
  int adapt_end_iter = (adaptation_end < 0) ?
    (adaptation_metric == "diversity" ? runs : burn) : adaptation_end;
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
  if (n<1||J<1||P<1||X.n_rows!=Y.n_rows||!X.is_finite()||Y.min()<0||
      s2.n_elem!=Y.n_cols||r2.n_elem!=Y.n_cols||!s2.is_finite()||
      !r2.is_finite()||s2.min()<=0||r2.min()<=0||a<=0||b<=0||c<=0)
    Rcpp::stop("Invalid DM data or prior parameters");
  for (int i=0;i<n;++i)
    if (arma::accu(arma::conv_to<arma::rowvec>::from(Y.row(i)))>
        std::numeric_limits<int>::max()) Rcpp::stop("DM row count exceeds integer range");
  dm_numeric_diagnostics::rejected_ratios=0;
  dm_numeric_diagnostics::precision_projections=0;
  dm_numeric_diagnostics::optimizer_fallbacks=0;

  if (target_effective_neighbors <= 0.0) {
    target_effective_neighbors = std::min(
      static_cast<double>(P), std::max(5.0, 0.025 * static_cast<double>(P))
    );
  }
  if (target_effective_neighbors < 1.0 ||
      target_effective_neighbors > static_cast<double>(P)) {
    Rcpp::stop("target_effective_neighbors must lie in [1, P].");
  }
  
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
  if (initial_xi.isNotNull()) {
    Rcpp::IntegerMatrix initial_xi_matrix(initial_xi);
    if (initial_xi_matrix.nrow() != P || initial_xi_matrix.ncol() != J) {
      Rcpp::stop("initial_xi must have dimensions P by J.");
    }
    xi = Rcpp::as<arma::imat>(initial_xi_matrix);
    if (arma::any(arma::vectorise(xi) < 0) ||
        arma::any(arma::vectorise(xi) > 1)) {
      Rcpp::stop("initial_xi must contain only zeros and ones.");
    }
  }
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
    log_prior_xi(j) = dm_indicator_log_prior(xi.col(j),a,b);
    lambda(j) = c * static_cast<double>(PP(j)) / static_cast<double>(n);
    if (PP(j) > 0) {
      beta_hat_memory[j][config[j]] = beta[j];
    }
  }
  
  arma::mat alpha_mat(n, J);
  alpha_mat.each_row() = alpha.t();
  arma::mat eta = alpha_mat + Xbeta;
  arma::mat gamma = arma::trunc_exp(eta);
  arma::vec gamma_sum_i = arma::sum(gamma, 1);
  
  double log_lik=dm_stable_log_kernel(Y,eta);
  if (!std::isfinite(log_lik)) Rcpp::stop("Nonfinite initial DM likelihood");
  double log_prior_alpha = -0.5 * arma::accu(arma::square(alpha) / s2);
  
  arma::mat prop_unif(P, J, arma::fill::ones);
  for (int j = 0; j < J; j++){
    prop_unif.col(j) /= P;
  }
  
  int save_every = (thin == 0) ? 1 : thin;
  int post_runs = (runs - burn + save_every - 1) / save_every;
  
  arma::mat alpha_post(post_runs, J, arma::fill::zeros);
  Rcpp::List beta_post(post_runs);
  arma::icube xi_post(P, J, post_runs, arma::fill::zeros);
  
  double s_alpha = 0.0;
  
  arma::vec m_alpha = alpha;
  
  arma::mat R_alpha = arma::eye(J, J); // 0.0005 *
  
  arma::mat Sigma_alpha = std::exp(s_alpha) * R_alpha;
  
  int accept_alpha = 0;
  arma::ivec joint_proposed(J,arma::fill::zeros),joint_accepted(J,arma::fill::zeros);
  arma::ivec within_proposed(J, arma::fill::zeros), within_accepted(J, arma::fill::zeros);
  arma::ivec within_fixed_proposed(J,arma::fill::zeros),within_fixed_accepted(J,arma::fill::zeros);
  arma::ivec within_curvature_proposed(J,arma::fill::zeros),within_curvature_accepted(J,arma::fill::zeros);
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
  arma::vec proposal_effective_window_sum(J, arma::fill::zeros);
  arma::vec proposal_effective_post_sum(J, arma::fill::zeros);
  arma::vec proposal_entropy_post_sum(J, arma::fill::zeros);
  arma::vec proposal_max_post_sum(J, arma::fill::zeros);
  arma::ivec proposal_post_count(J, arma::fill::zeros);
  arma::vec proposal_effective_this_iter(J, arma::fill::zeros);
  arma::vec proposal_entropy_this_iter(J, arma::fill::zeros);
  arma::vec proposal_max_this_iter(J, arma::fill::zeros);
  int max_adapt_epochs = 0;
  if (adapt_g) {
    max_adapt_epochs = (adapt_end_iter - adapt_start_iter + 1) / window_size;
  }
  arma::mat g_window_acceptance(max_adapt_epochs, J);
  arma::mat g_window_values(max_adapt_epochs, J);
  arma::mat g_window_effective_neighbors(max_adapt_epochs, J);
  g_window_acceptance.fill(arma::datum::nan);
  g_window_values.fill(arma::datum::nan);
  g_window_effective_neighbors.fill(arma::datum::nan);

  arma::mat g_trace(post_runs, J, arma::fill::value(arma::datum::nan));
  arma::mat proposal_effective_trace(post_runs, J, arma::fill::value(arma::datum::nan));
  arma::mat proposal_entropy_trace(post_runs, J, arma::fill::value(arma::datum::nan));
  arma::mat proposal_max_trace(post_runs, J, arma::fill::value(arma::datum::nan));
  arma::vec log_posterior_trace(post_runs, arma::fill::value(arma::datum::nan));
  
  double log_g_min = std::log(min_g);
  double log_g_max = std::log(max_g);
  for (int j = 0; j < J; j++) {
    g(j) = clamp_double(g_initial, min_g, max_g);
    log_g(j) = std::log(g(j));
    prev_log_g(j) = log_g(j);
    prev_prev_log_g(j) = log_g(j);
  }
  
  // Respect the OpenMP runtime configuration.  On Hoffman2 the job scripts set
  // OMP_NUM_THREADS to the Grid Engine allocation (NSLOTS); a hard-coded value
  // here would oversubscribe the node and invalidate timing comparisons.
  
// Console-safety-only regression build: never call R output on an OpenMP worker.
int dm_banner_threads = 1;
#pragma omp parallel
{
#pragma omp master
{
  dm_banner_threads = omp_get_num_threads();
}
}
Rcpp::Rcout << "### Using " << dm_banner_threads << " threads \n" << std::endl;

// Segment boundaries are whole sweeps, never new warmup schedules.
const int stop_iter=stop_after<0 ? runs : stop_after;
if(stop_iter<1 || stop_iter>runs)Rcpp::stop("Invalid checkpoint endpoint");
int completed_before=0;
double elapsed_before=0;
Rcpp::List checkpoint_config;
checkpoint_config["runs"]=runs;
checkpoint_config["burn"]=burn;
checkpoint_config["thin"]=thin;
checkpoint_config["informed"]=informed;
checkpoint_config["trgt_acpt_alpha"]=trgt_acpt_alpha;
checkpoint_config["adapt_g"]=adapt_g;
checkpoint_config["adaptation_start"]=adaptation_start;
checkpoint_config["adaptation_end"]=adaptation_end;
checkpoint_config["window_size"]=window_size;
checkpoint_config["rm_scale"]=rm_scale;
checkpoint_config["decay_power"]=decay_power;
checkpoint_config["min_g"]=min_g;
checkpoint_config["max_g"]=max_g;
checkpoint_config["use_move"]=use_move;
checkpoint_config["n_move_attempts"]=n_move_attempts;
checkpoint_config["lambda_move"]=lambda_move;
checkpoint_config["adaptation_metric"]=adaptation_metric;
checkpoint_config["target_effective_neighbors"]=target_effective_neighbors;
checkpoint_config["adaptation_offset"]=adaptation_offset;
checkpoint_config["balance_directions"]=balance_directions;
checkpoint_config["uniform_mixture"]=uniform_mixture;
checkpoint_config["score_optimizer_maxit"]=score_optimizer_maxit;
checkpoint_config["coefficient_optimizer_maxit"]=coefficient_optimizer_maxit;
checkpoint_config["cache_score_fits"]=cache_score_fits;
checkpoint_config["discrete_proposal"]=discrete_proposal;
checkpoint_config["laplace_optimizer_maxit"]=laplace_optimizer_maxit;
checkpoint_config["coefficient_warm_start"]=coefficient_warm_start;
checkpoint_config["within_model_update"]=within_model_update;
checkpoint_config["within_model_attempts"]=within_model_attempts;
checkpoint_config["within_model_scale"]=within_model_scale;
checkpoint_config["within_model_fixed_probability"]=within_model_fixed_probability;
checkpoint_config["within_model_fixed_scale"]=within_model_fixed_scale;
checkpoint_config["X"]=X;
checkpoint_config["Y"]=Y;
checkpoint_config["s2"]=s2;
checkpoint_config["r2"]=r2;
checkpoint_config["a"]=a;
checkpoint_config["b"]=b;
checkpoint_config["g_initial"]=g_initial;
checkpoint_config["c"]=c;
checkpoint_config["neigh"]=neigh.isNotNull() ? Rcpp::List(neigh) : Rcpp::List();
checkpoint_config["initial_xi"]=initial_xi.isNotNull() ? Rcpp::IntegerMatrix(initial_xi) : Rcpp::IntegerMatrix();
checkpoint_config["threads"]=omp_get_max_threads();
if(resume_state.isNotNull()) {
  if(resume_draws.isNull())Rcpp::stop("Resume needs retained draws AND full continuation state");
  Rcpp::List state(resume_state),draws(resume_draws);
  Rcpp::Function identical("identical");
  if(Rcpp::as<int>(state["schema"])!=1 ||
     !Rcpp::as<bool>(identical(state["config"],checkpoint_config)))
    Rcpp::stop("Checkpoint configuration/data/thread count differs");
  completed_before=Rcpp::as<int>(state["completed"]);
  elapsed_before=Rcpp::as<double>(state["elapsed_seconds"]);
  if(completed_before<1 || completed_before>=stop_iter || !std::isfinite(elapsed_before))
    Rcpp::stop("Invalid saved checkpoint boundary");
  alpha=Rcpp::as<arma::vec>(state["alpha"]);
  log_prior_beta=Rcpp::as<arma::vec>(state["log_prior_beta"]);
  log_prior_xi=Rcpp::as<arma::vec>(state["log_prior_xi"]);
  lambda=Rcpp::as<arma::vec>(state["lambda"]);
  gamma_sum_i=Rcpp::as<arma::vec>(state["gamma_sum_i"]);
  m_alpha=Rcpp::as<arma::vec>(state["m_alpha"]);
  log_g=Rcpp::as<arma::vec>(state["log_g"]);
  prev_log_g=Rcpp::as<arma::vec>(state["prev_log_g"]);
  prev_prev_log_g=Rcpp::as<arma::vec>(state["prev_prev_log_g"]);
  g=Rcpp::as<arma::vec>(state["g"]);
  prev_accept_rate_g=Rcpp::as<arma::vec>(state["prev_accept_rate_g"]);
  proposal_effective_window_sum=Rcpp::as<arma::vec>(state["proposal_effective_window_sum"]);
  proposal_effective_post_sum=Rcpp::as<arma::vec>(state["proposal_effective_post_sum"]);
  proposal_entropy_post_sum=Rcpp::as<arma::vec>(state["proposal_entropy_post_sum"]);
  proposal_max_post_sum=Rcpp::as<arma::vec>(state["proposal_max_post_sum"]);
  proposal_effective_this_iter=Rcpp::as<arma::vec>(state["proposal_effective_this_iter"]);
  proposal_entropy_this_iter=Rcpp::as<arma::vec>(state["proposal_entropy_this_iter"]);
  proposal_max_this_iter=Rcpp::as<arma::vec>(state["proposal_max_this_iter"]);
  xi=Rcpp::as<arma::imat>(state["xi"]);
  PP=Rcpp::as<arma::ivec>(state["PP"]);
  within_proposed=Rcpp::as<arma::ivec>(state["within_proposed"]);
  within_accepted=Rcpp::as<arma::ivec>(state["within_accepted"]);
  within_fixed_proposed=Rcpp::as<arma::ivec>(state["within_fixed_proposed"]);
  within_fixed_accepted=Rcpp::as<arma::ivec>(state["within_fixed_accepted"]);
  within_curvature_proposed=Rcpp::as<arma::ivec>(state["within_curvature_proposed"]);
  within_curvature_accepted=Rcpp::as<arma::ivec>(state["within_curvature_accepted"]);
  accept_betaxi=Rcpp::as<arma::ivec>(state["accept_betaxi"]);
  accept_move=Rcpp::as<arma::ivec>(state["accept_move"]);
  prop_move=Rcpp::as<arma::ivec>(state["prop_move"]);
  acc_betaxi_window=Rcpp::as<arma::ivec>(state["acc_betaxi_window"]);
  iter_counter_g=Rcpp::as<arma::ivec>(state["iter_counter_g"]);
  update_counter_g=Rcpp::as<arma::ivec>(state["update_counter_g"]);
  proposal_post_count=Rcpp::as<arma::ivec>(state["proposal_post_count"]);
  Xbeta=Rcpp::as<arma::mat>(state["Xbeta"]);
  alpha_mat=Rcpp::as<arma::mat>(state["alpha_mat"]);
  eta=Rcpp::as<arma::mat>(state["eta"]);
  gamma=Rcpp::as<arma::mat>(state["gamma"]);
  R_alpha=Rcpp::as<arma::mat>(state["R_alpha"]);
  Sigma_alpha=Rcpp::as<arma::mat>(state["Sigma_alpha"]);
  g_window_acceptance=Rcpp::as<arma::mat>(state["g_window_acceptance"]);
  g_window_values=Rcpp::as<arma::mat>(state["g_window_values"]);
  g_window_effective_neighbors=Rcpp::as<arma::mat>(state["g_window_effective_neighbors"]);
  log_lik=Rcpp::as<double>(state["log_lik"]);
  log_prior_alpha=Rcpp::as<double>(state["log_prior_alpha"]);
  s_alpha=Rcpp::as<double>(state["s_alpha"]);
  accept_alpha=Rcpp::as<int>(state["accept_alpha"]);
  accept_betaxi_global=Rcpp::as<int>(state["accept_betaxi_global"]);
  accept_move_global=Rcpp::as<int>(state["accept_move_global"]);
  prop_move_global=Rcpp::as<int>(state["prop_move_global"]);
  beta=Rcpp::as<std::vector<arma::vec>>(state["beta"]);
  dm_numeric_diagnostics::rejected_ratios=Rcpp::as<int>(state["rejected_ratios"]);
  dm_numeric_diagnostics::precision_projections=Rcpp::as<int>(state["precision_projections"]);
  dm_numeric_diagnostics::optimizer_fallbacks=Rcpp::as<int>(state["optimizer_fallbacks"]);
  if(alpha.n_elem!=static_cast<arma::uword>(J)||xi.n_rows!=static_cast<arma::uword>(P)||
     xi.n_cols!=static_cast<arma::uword>(J)||beta.size()!=static_cast<size_t>(J)||
     Xbeta.n_rows!=static_cast<arma::uword>(n)||Xbeta.n_cols!=static_cast<arma::uword>(J))
    Rcpp::stop("Checkpoint latent dimensions differ");
  // Rebuild deterministic views only; restore Xbeta/eta/likelihood exactly, not
  // by new multiplication or summation. Historical fit caches affect cost only.
  for(int j=0;j<J;++j) {
    config[j]=vectorToKey(arma::conv_to<arma::ivec>::from(xi.col(j)));
    selected[j]=arma::find(xi.col(j)>0);X_list[j]=X.cols(selected[j]);
    if(PP(j)!=static_cast<int>(selected[j].n_elem)||beta[j].n_elem!=selected[j].n_elem)
      Rcpp::stop("Checkpoint coefficient/model dimensions differ");
    beta_hat_memory[j].clear();
    if(PP(j)>0)beta_hat_memory[j][config[j]]=beta[j];
  }
  const int kept=completed_before<=burn ? 0 : (completed_before-burn+save_every-1)/save_every;
  arma::mat old_alpha=Rcpp::as<arma::mat>(draws["alpha"]);
  arma::icube old_xi=Rcpp::as<arma::icube>(draws["xi"]);
  Rcpp::List old_beta=draws["beta"];
  if(old_alpha.n_rows!=static_cast<arma::uword>(kept)||old_alpha.n_cols!=static_cast<arma::uword>(J)||
     old_xi.n_rows!=static_cast<arma::uword>(P)||old_xi.n_cols!=static_cast<arma::uword>(J)||
     old_xi.n_slices!=static_cast<arma::uword>(kept)||old_beta.size()!=kept)
    Rcpp::stop("Checkpoint retained dimensions differ");
  if(kept>0) {
    alpha_post.rows(0,kept-1)=old_alpha;xi_post.slices(0,kept-1)=old_xi;
    for(int d=0;d<kept;++d)beta_post[d]=old_beta[d];
    g_trace.rows(0,kept-1)=Rcpp::as<arma::mat>(draws["g_trace"]);
    proposal_effective_trace.rows(0,kept-1)=Rcpp::as<arma::mat>(draws["proposal_effective_neighbors_trace"]);
    proposal_entropy_trace.rows(0,kept-1)=Rcpp::as<arma::mat>(draws["proposal_entropy_trace"]);
    proposal_max_trace.rows(0,kept-1)=Rcpp::as<arma::mat>(draws["proposal_max_probability_trace"]);
    log_posterior_trace.head(kept)=Rcpp::as<arma::vec>(draws["log_posterior_trace"]);
  }
} else if(resume_draws.isNotNull())Rcpp::stop("Retained draws without continuation state");

// A NEW diagnostic trajectory, not a resumed/extended production chain.
// No RNG is consumed by this initialization. Only target quantities and views
// are recomputed; the old trained intercept covariance and g remain fixed.
if(pilot_start.isNotNull()) {
  Rcpp::List start(pilot_start);
  alpha=Rcpp::as<arma::vec>(start["alpha"]);xi=Rcpp::as<arma::imat>(start["xi"]);
  beta=Rcpp::as<std::vector<arma::vec>>(start["beta"]);
  g=Rcpp::as<arma::vec>(start["g"]);Sigma_alpha=Rcpp::as<arma::mat>(start["Sigma_alpha"]);
  if(alpha.n_elem!=static_cast<arma::uword>(J)||xi.n_rows!=static_cast<arma::uword>(P)||
      xi.n_cols!=static_cast<arma::uword>(J)||beta.size()!=static_cast<size_t>(J)||
      g.n_elem!=static_cast<arma::uword>(J)||!alpha.is_finite()||!g.is_finite()||
      g.min()<=0||Sigma_alpha.n_rows!=static_cast<arma::uword>(J)||
      Sigma_alpha.n_cols!=static_cast<arma::uword>(J)||!Sigma_alpha.is_finite()||
      arma::any(arma::vectorise(xi)<0)||arma::any(arma::vectorise(xi)>1))
    Rcpp::stop("Invalid saved-endpoint pilot initialization");
  arma::mat check_chol;
  if(!arma::chol(check_chol,arma::symmatu(Sigma_alpha)))Rcpp::stop("Pilot alpha covariance is not SPD");
  PP=arma::sum(xi,0).t();
  for(int j=0;j<J;++j) {
    selected[j]=arma::find(xi.col(j)>0);X_list[j]=X.cols(selected[j]);
    config[j]=vectorToKey(arma::conv_to<arma::ivec>::from(xi.col(j)));
    if(beta[j].n_elem!=selected[j].n_elem||!beta[j].is_finite())Rcpp::stop("Pilot beta/indicator mismatch");
    Xbeta.col(j)=X_list[j]*beta[j];
    log_prior_beta(j)=-(PP(j)/2.0)*(std::log(2.0*M_PI)+std::log(r2(j)))-.5*arma::dot(beta[j],beta[j])/r2(j);
    log_prior_xi(j)=dm_indicator_log_prior(xi.col(j),a,b);
    lambda(j)=c*static_cast<double>(PP(j))/static_cast<double>(n);
    beta_hat_memory[j].clear();if(PP(j)>0)beta_hat_memory[j][config[j]]=beta[j];
    log_g(j)=std::log(g(j));prev_log_g(j)=log_g(j);prev_prev_log_g(j)=log_g(j);
  }
  alpha_mat.each_row()=alpha.t();eta=alpha_mat+Xbeta;gamma=arma::trunc_exp(eta);
  gamma_sum_i=arma::sum(gamma,1);log_lik=dm_stable_log_kernel(Y,eta);
  log_prior_alpha=-.5*arma::accu(arma::square(alpha)/s2);
}
const double pilot_initial_log_posterior=log_lik+log_prior_alpha+
  arma::accu(log_prior_beta)+arma::accu(log_prior_xi);

auto chain_start = std::chrono::steady_clock::now();
int next_progress_hundredth_percent = completed_before==0 ? 1 :
  ((static_cast<int>(10000.0*completed_before/runs)/500)+1)*500;

for (int iter = completed_before; iter < stop_iter; iter++) {
  if(iter%25==0)Rcpp::checkUserInterrupt();
  int iter_id = iter + 1;
  if (!std::isfinite(log_lik)||!std::isfinite(log_prior_alpha)||
      !log_prior_beta.is_finite()||!log_prior_xi.is_finite()||
      !alpha.is_finite()||!Xbeta.is_finite())
    Rcpp::stop("Nonfinite current DM state at iteration %d",iter_id);
  
  double iter07 = pow(iter + 2, -0.7);
  
  //////////////////
  // alpha update //
  //////////////////
  
  if (!Sigma_alpha.is_finite()) Rcpp::stop("Nonfinite alpha covariance at iteration %d",iter_id);
  Sigma_alpha = arma::symmatu(Sigma_alpha);
  arma::vec alpha_prop = arma::mvnrnd(alpha, Sigma_alpha);
  arma::mat alpha_mat_prop(n, J);
  alpha_mat_prop.each_row() = alpha_prop.t();
  
  arma::mat eta_prop = alpha_mat_prop + Xbeta;
  arma::mat gamma_prop = arma::trunc_exp(eta_prop);
  arma::vec gamma_prop_sum_i = arma::sum(gamma_prop, 1);
  
  double log_lik_prop=dm_stable_log_kernel(Y,eta_prop);

  double log_prior_alpha_prop = -0.5 * arma::accu(arma::square(alpha_prop) / s2);
  
  double log_ratio_alpha = log_lik_prop - log_lik + log_prior_alpha_prop - log_prior_alpha;
  
  double u = arma::randu();
  
  double log_acpt_alpha = dm_safe_log_acceptance(log_ratio_alpha);
  
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
  
  // Freeze the intercept proposal as well as g for retained sampling.
  if (iter < burn) {
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
  Sigma_alpha = arma::symmatu(Sigma_alpha);
  }
  
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
      if (discrete_proposal == "gagnon_laplace") {
        proposal_xi_j = proposal_gagnon_laplace(
          xi.col(j), j, PP, Xbeta, beta_hat_memory[j], Y, X, alpha,
          alpha_mat, y_sum_i, log_lik, gamma, a, b, r2(j),
          laplace_optimizer_maxit, cache_score_fits
        );
      } else {
        proposal_xi_j = proposal_inf(
          xi.col(j), j, PP, Xbeta, beta_hat_memory[j], Y, X, alpha,
          alpha_mat, y_sum_i, log_lik, gamma, g(j), c,
          balance_directions, uniform_mixture, score_optimizer_maxit,
          cache_score_fits
        );
      }
    } else {
      proposal_xi_j = prop_unif.col(j);  
    }

    double proposal_square_sum = arma::accu(arma::square(proposal_xi_j));
    proposal_effective_this_iter(j) = 1.0 / proposal_square_sum;
    arma::uvec proposal_positive = arma::find(proposal_xi_j > 0.0);
    proposal_entropy_this_iter(j) = -arma::accu(
      proposal_xi_j.elem(proposal_positive) %
        arma::log(proposal_xi_j.elem(proposal_positive))
    );
    proposal_max_this_iter(j) = proposal_xi_j.max();
    
    int ind_j = sample_int(P, 1, proposal_xi_j) - 1;
    xi_prop.col(j) = xi.col(j);
    xi_prop(ind_j, j) = 1 - xi(ind_j, j);
    
    config_prop[j] = vectorToKey(arma::conv_to<arma::ivec>::from(xi_prop.col(j)));
    
    PP_prop = PP;
    PP_prop(j) = arma::sum(xi_prop.col(j));
    
    selected_prop[j] = arma::find(xi_prop.col(j) > 0);
    
    X_list_prop[j] = X.cols(selected_prop[j]);
    
    log_prior_xi_prop(j) = dm_indicator_log_prior(xi_prop.col(j),a,b);
    
    lambda_prop(j) = c * static_cast<double>(PP_prop(j)) / static_cast<double>(n);
    
    // Proposed state proposal (only if PP_prop(j) > 0)
    arma::vec beta_j_prop_hat;
    arma::mat H_beta_j_prop;
    double log_proposal_beta_forward;
    
    Xbeta_prop = Xbeta;
    
    if (PP_prop(j) > 0) {
      auto it = beta_hat_memory[j].find(config_prop[j]);
      double proposal_penalty = discrete_proposal == "gagnon_laplace" ?
        1.0 / r2(j) : lambda_prop(j);
      if (coefficient_warm_start) {
        arma::vec start=dm_project_state_start(selected[j],beta[j],selected_prop[j]);
        beta_j_prop_hat=dm_reg_lbfgs(X_list_prop[j],Y,alpha.t(),PP_prop,Xbeta,j,
          proposal_penalty,start,coefficient_optimizer_maxit,true);
      } else if (it != beta_hat_memory[j].end()) {
        // Warm start from previous beta_hat
        beta_j_prop_hat = dm_reg_lbfgs(X_list_prop[j], Y, alpha.t(), PP_prop, Xbeta, j, proposal_penalty, it->second, coefficient_optimizer_maxit);
      } else {
        // New configuration: start from random
        arma::vec beta_start(PP_prop(j), arma::fill::zeros);
        beta_j_prop_hat = dm_reg_lbfgs(X_list_prop[j], Y, alpha.t(), PP_prop, Xbeta, j, proposal_penalty, beta_start, coefficient_optimizer_maxit);
      }
      // Store optimized beta_hat in memory
      beta_hat_memory[j][config_prop[j]] = beta_j_prop_hat;
      
      arma::mat H_beta_j_prop = dm_reg_hessian(X_list_prop[j], Y, alpha.t(), PP_prop, Xbeta, j, beta_j_prop_hat, proposal_penalty);
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
      double reverse_proposal_penalty = discrete_proposal == "gagnon_laplace" ?
        1.0 / r2(j) : lambda(j);
      arma::vec reverse_start=coefficient_warm_start ?
        dm_project_state_start(selected_prop[j],beta_prop[j],selected[j]) : beta_hat_memory[j][config[j]];
      arma::vec beta_j_hat_reverse = dm_reg_lbfgs(X_list[j], Y, alpha.t(), PP, Xbeta_prop, j, reverse_proposal_penalty, reverse_start, coefficient_optimizer_maxit, coefficient_warm_start);
      arma::mat H_beta_j_reverse = dm_reg_hessian(X_list[j], Y, alpha.t(), PP, Xbeta_prop, j, beta_j_hat_reverse, reverse_proposal_penalty);
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
    
    double log_lik_prop=log_lik+
      dm_stable_replacement_delta(Y,alpha_mat,Xbeta,j,Xbeta_prop.col(j));

    u = arma::randu();
    
    arma::vec proposal_xi_j_prop(P);
    
    if (informed == true) {
      if (discrete_proposal == "gagnon_laplace") {
        proposal_xi_j_prop = proposal_gagnon_laplace(
          xi_prop.col(j), j, PP_prop, Xbeta_prop, beta_hat_memory[j], Y, X,
          alpha, alpha_mat, y_sum_i, log_lik_prop, gamma_prop, a, b,
          r2(j), laplace_optimizer_maxit, cache_score_fits
        );
      } else {
        proposal_xi_j_prop = proposal_inf(
          xi_prop.col(j), j, PP_prop, Xbeta_prop, beta_hat_memory[j], Y, X,
          alpha, alpha_mat, y_sum_i, log_lik_prop, gamma_prop, g(j), c,
          balance_directions, uniform_mixture, score_optimizer_maxit,
          cache_score_fits
        );
      }
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
    
    if (log(u) < dm_safe_log_acceptance(log_ratio_betaxi_j)) {
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
          log_lik, gamma, lambda_move, c, score_optimizer_maxit
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
        
        double log_prior_xi_move = dm_indicator_log_prior(xi_move_j,a,b);
        double lambda_move_j =
          c * static_cast<double>(PP_move(j)) / static_cast<double>(n);
        
        arma::vec beta_j_move_hat;
        double log_proposal_beta_forward = 0.0;
        arma::vec beta_move_j;
        arma::mat Xbeta_move = Xbeta;
        
        if (PP_move(j) > 0) {
          auto it_move = beta_hat_memory[j].find(config_move);
          if (coefficient_warm_start) {
            arma::vec start=dm_project_state_start(selected[j],beta[j],selected_move);
            beta_j_move_hat=dm_reg_lbfgs(X_move,Y,alpha.t(),PP_move,Xbeta,j,
              lambda_move_j,start,coefficient_optimizer_maxit,true);
          } else if (it_move != beta_hat_memory[j].end()) {
            beta_j_move_hat = dm_reg_lbfgs(
              X_move, Y, alpha.t(), PP_move, Xbeta, j, lambda_move_j,
              it_move->second, coefficient_optimizer_maxit
            );
          } else {
            arma::vec beta_start(PP_move(j), arma::fill::zeros);
            beta_j_move_hat = dm_reg_lbfgs(
              X_move, Y, alpha.t(), PP_move, Xbeta, j, lambda_move_j,
              beta_start, coefficient_optimizer_maxit
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
        
        double log_lik_move=log_lik+
          dm_stable_replacement_delta(Y,alpha_mat,Xbeta,j,Xbeta_move.col(j));

        double log_proposal_beta_reverse = 0.0;
        if (PP(j) > 0) {
          auto it_reverse = beta_hat_memory[j].find(config[j]);
          arma::vec beta_start_reverse = (it_reverse != beta_hat_memory[j].end()) ?
          it_reverse->second : beta[j];
          if (coefficient_warm_start)
            beta_start_reverse=dm_project_state_start(selected_move,beta_move_j,selected[j]);
          arma::vec beta_j_hat_reverse = dm_reg_lbfgs(
            X_list[j], Y, alpha.t(), PP, Xbeta_move, j, lambda(j),
            beta_start_reverse, coefficient_optimizer_maxit, coefficient_warm_start
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
          log_lik_move, gamma_move, lambda_move, c, score_optimizer_maxit
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
        
        if (std::log(arma::randu()) < dm_safe_log_acceptance(log_ratio_move)) {
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
    
    // Optional fixed-weight mixture of within-model RW kernels. No PMLE fitting.
    // Curvature recomputes reverse geometry; design-only covariance is symmetric.
    if (within_model_update && PP(j)>0) {
      for (int attempt=0;attempt<within_model_attempts;++attempt) {
        ++within_proposed(j);
        // Fixed state-independent selection gives a mixture of two MH kernels.
        // Probability zero consumes NO extra random variate (legacy equivalence).
        bool fixed_design=within_model_fixed_probability==1.0 ||
          (within_model_fixed_probability>0 && arma::randu()<within_model_fixed_probability);
        if (fixed_design) ++within_fixed_proposed(j); else ++within_curvature_proposed(j);
        double step_scale=fixed_design ? within_model_fixed_scale : within_model_scale;
        arma::mat raw= fixed_design ? arma::mat(X_list[j].t()*X_list[j]/Y.n_rows+
          arma::eye(PP(j),PP(j))/r2(j)) : arma::mat(-dm_reg_hessian(
          X_list[j],Y,alpha.t(),PP,Xbeta,j,beta[j],1.0/r2(j)));
        arma::mat qf=dm_gaussian_precision(raw)/(step_scale*step_scale);
        arma::vec candidate=rmvn_prec(beta[j],qf);
        arma::mat candidate_xb=Xbeta;
        candidate_xb.col(j)=X_list[j]*candidate;
        DMWithinTerms terms=dm_within_terms(X_list[j],Y,alpha,PP,Xbeta,j,beta[j],
          candidate,r2(j),step_scale,qf,fixed_design);
        if (std::log(arma::randu())<dm_safe_log_acceptance(terms.ratio)) {
          beta[j]=candidate;
          Xbeta.col(j)=candidate_xb.col(j);
          eta.col(j)=alpha_mat.col(j)+Xbeta.col(j);
          arma::vec new_gamma=arma::trunc_exp(eta.col(j));
          gamma_sum_i+=new_gamma-gamma.col(j);
          gamma.col(j)=new_gamma;
          log_lik+=terms.likelihood_delta;
          log_prior_beta(j)+=terms.prior_delta;
          ++within_accepted(j);
          if (fixed_design) ++within_fixed_accepted(j); else ++within_curvature_accepted(j);
        }
      }
    }

    // Optional ADDITIONAL full-rank joint MH kernel, shared by all methods.
    // All existing updates stay unchanged. Disabled consumes NO random draws.
    if(joint_update) for(int attempt=0;attempt<joint_attempts;++attempt) {
      ++joint_proposed(j);
      arma::vec current=dm_joint::coordinates(alpha(j),beta[j]);
      arma::mat raw=dm_joint::raw_precision(X_list[j],Y,alpha,PP,Xbeta,j,current,s2(j),r2(j));
      arma::mat qf=dm_joint::precision(raw,joint_scale);
      arma::vec candidate=rmvn_prec(current,qf);
      dm_joint::Terms terms=dm_joint::terms(X_list[j],Y,alpha,PP,Xbeta,j,beta[j],candidate,
        s2(j),r2(j),joint_scale,raw,qf);
      if(std::log(arma::randu())<dm_safe_log_acceptance(terms.ratio)) {
        alpha(j)=candidate(0);beta[j]=candidate.tail(PP(j));
        alpha_mat.col(j).fill(alpha(j));Xbeta.col(j)=X_list[j]*beta[j];
        eta.col(j)=alpha_mat.col(j)+Xbeta.col(j);
        arma::vec new_gamma=arma::trunc_exp(eta.col(j));
        gamma_sum_i+=new_gamma-gamma.col(j);gamma.col(j)=new_gamma;
        log_lik+=terms.likelihood_delta;log_prior_alpha+=terms.alpha_prior_delta;
        log_prior_beta(j)+=terms.beta_prior_delta;++joint_accepted(j);
        // All cached conditioning keys include alpha and other-category Xbeta;
        // the next score/fit therefore cannot reuse an obsolete conditioning state.
      }
    }

    bool in_adaptation_window = adapt_g &&
      (iter_id >= adapt_start_iter) &&
      (iter_id <= adapt_end_iter);
    if (in_adaptation_window) {
      acc_betaxi_window(j) += (accepted_betaxi_j ? 1 : 0);
      proposal_effective_window_sum(j) += proposal_effective_this_iter(j);
      iter_counter_g(j)++;
    }

    if (iter >= burn) {
      proposal_effective_post_sum(j) += proposal_effective_this_iter(j);
      proposal_entropy_post_sum(j) += proposal_entropy_this_iter(j);
      proposal_max_post_sum(j) += proposal_max_this_iter(j);
      proposal_post_count(j)++;
    }
    
  }
  
  // After the loop: increment counter if any beta changed
  if (any_betaxi_changed) {
    accept_betaxi_global++;  // Increment a single counter (not indexed by j)
  }
  
  ////////////////////////
  // DIMINISHING ONLINE ADAPTATION FOR g (category-wise)
  ////////////////////////
  
  // This should be placed after the for loop over j (after all beta, xi updates)
  for (int j = 0; j < J; j++) {
    
    if (adapt_g && iter_counter_g(j) == window_size) {
      update_counter_g(j) += 1;
      double curr_accept_rate_g = static_cast<double>(acc_betaxi_window(j)) /
        static_cast<double>(window_size);
      double step_size = rm_scale *
        std::pow(static_cast<double>(update_counter_g(j)), -decay_power);
      double mean_effective_neighbors = proposal_effective_window_sum(j) /
        static_cast<double>(window_size);
      
      // Bootstrap one bounded move after the first window, then use
      // the Robbins-Monro hill-climbing update thereafter.
      if (adaptation_metric == "diversity") {
        double gain_index = std::max(
          1.0, adaptation_offset + static_cast<double>(update_counter_g(j) - 1)
        );
        double diversity_step = rm_scale * std::pow(gain_index, -decay_power);
        log_g(j) += diversity_step * (
          std::log(mean_effective_neighbors) -
            std::log(target_effective_neighbors)
        );
        log_g(j) = clamp_double(log_g(j), log_g_min, log_g_max);
        g(j) = std::exp(log_g(j));
      } else if (update_counter_g(j) == 1) {
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
        g_window_effective_neighbors(epoch_idx, j) = mean_effective_neighbors;
      }
      prev_accept_rate_g(j) = curr_accept_rate_g;
      acc_betaxi_window(j) = 0;
      proposal_effective_window_sum(j) = 0.0;
      iter_counter_g(j) = 0;
    }
  }
  
  // In the loop
  if((iter >= burn) && ((iter - burn) % save_every == 0)){
    int save_idx = (iter - burn) / save_every;
    alpha_post.row(save_idx) = alpha.t();
    beta_post[save_idx] = Rcpp::wrap(beta);  // Wrap converts std::vector<arma::vec> to R list
    xi_post.slice(save_idx) = xi;
    g_trace.row(save_idx) = g.t();
    proposal_effective_trace.row(save_idx) = proposal_effective_this_iter.t();
    proposal_entropy_trace.row(save_idx) = proposal_entropy_this_iter.t();
    proposal_max_trace.row(save_idx) = proposal_max_this_iter.t();
    // Refresh from the complete state to prevent accumulated delta drift.
    log_lik=dm_stable_log_kernel(Y,alpha_mat+Xbeta);
    if (!std::isfinite(log_lik)) Rcpp::stop("Nonfinite retained DM likelihood");
    log_posterior_trace(save_idx) = log_lik + log_prior_alpha +
      arma::accu(log_prior_beta) + arma::accu(log_prior_xi);
    if (!std::isfinite(log_posterior_trace(save_idx)))
      Rcpp::stop("Nonfinite retained DM posterior at iteration %d",iter_id);
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
      elapsed_before+std::chrono::duration<double>(std::chrono::steady_clock::now() - chain_start).count();
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
    // R console callbacks / PowerShell-captured stdout may discard carriage
    // returns. Write compact terminal updates directly to stderr instead.
#ifdef _WIN32
    const bool terminal_progress = _isatty(_fileno(stderr)) != 0;
#else
    const bool terminal_progress = isatty(fileno(stderr)) != 0;
#endif
    if (terminal_progress) {
      std::ostringstream compact;
      compact << std::fixed << std::setprecision(1) << progress << "%  "
        << elapsed_seconds << "s  ETA " << remaining_seconds << "s"
        << std::setprecision(2) << " | alpha " << acc_rate_alpha
        << " flip " << acc_rate_betaxi_global << " swap " << acc_rate_move_global;
      std::string line=compact.str();
      if (line.size()>75) line.resize(75);
      std::fprintf(stderr,"\r%-75s",line.c_str());
      std::fflush(stderr);
      if (iter+1==stop_iter) std::fprintf(stderr,"\n");
    } else {
      // Batch/captured output must have real line endings, not concatenated CRs.
      Rprintf("%s\n",progress_line.c_str());
      R_FlushConsole();
    }
    next_progress_hundredth_percent =
      ((completed_progress_hundredth_percent / 500) + 1) * 500;
  }
  
}

Rprintf("\n");
R_FlushConsole();

// Trim ONLY retained output buffers, not the full planned adaptation schedule.
const int kept=stop_iter<=burn ? 0 : (stop_iter-burn+save_every-1)/save_every;
alpha_post.resize(kept,J);xi_post.resize(P,J,kept);
g_trace.resize(kept,J);proposal_effective_trace.resize(kept,J);
proposal_entropy_trace.resize(kept,J);proposal_max_trace.resize(kept,J);
log_posterior_trace.resize(kept);
Rcpp::List beta_kept(kept);
for(int d=0;d<kept;++d)beta_kept[d]=beta_post[d];
beta_post=beta_kept;
arma::vec acpt_rate_betaxi_by_category =
  arma::conv_to<arma::vec>::from(accept_betaxi) / static_cast<double>(stop_iter);
arma::vec acpt_rate_move_by_category(J, arma::fill::zeros);
arma::vec mean_proposal_effective_neighbors(J, arma::fill::value(arma::datum::nan));
arma::vec mean_proposal_entropy(J, arma::fill::value(arma::datum::nan));
arma::vec mean_max_proposal_probability(J, arma::fill::value(arma::datum::nan));
for (int j = 0; j < J; ++j) {
  if (prop_move(j) > 0) {
    acpt_rate_move_by_category(j) =
      static_cast<double>(accept_move(j)) / static_cast<double>(prop_move(j));
  }
  if (proposal_post_count(j) > 0) {
    double denominator = static_cast<double>(proposal_post_count(j));
    mean_proposal_effective_neighbors(j) = proposal_effective_post_sum(j) / denominator;
    mean_proposal_entropy(j) = proposal_entropy_post_sum(j) / denominator;
    mean_max_proposal_probability(j) = proposal_max_post_sum(j) / denominator;
  }
}

double total_seconds = elapsed_before+std::chrono::duration<double>(
  std::chrono::steady_clock::now() - chain_start
).count();

Rcpp::List output = Rcpp::List::create(
  Rcpp::Named("alpha") = alpha_post,
  Rcpp::Named("beta") = beta_post,
  Rcpp::Named("xi") = xi_post,
  Rcpp::Named("acpt_rate_alpha") = static_cast<double>(accept_alpha) / static_cast<double>(stop_iter),
  Rcpp::Named("acpt_rate_betaxi_global") = static_cast<double>(accept_betaxi_global) / static_cast<double>(stop_iter),
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

output["g_window_effective_neighbors"] = g_window_effective_neighbors;
output["g_trace"] = g_trace;
output["proposal_effective_neighbors_trace"] = proposal_effective_trace;
output["proposal_entropy_trace"] = proposal_entropy_trace;
output["proposal_max_probability_trace"] = proposal_max_trace;
output["log_posterior_trace"] = log_posterior_trace;
output["mean_proposal_effective_neighbors"] = mean_proposal_effective_neighbors;
output["mean_proposal_entropy"] = mean_proposal_entropy;
output["mean_max_proposal_probability"] = mean_max_proposal_probability;
output["adaptation_updates_by_category"] = update_counter_g;
output["adaptation_metric"] = adaptation_metric;
output["target_effective_neighbors"] = target_effective_neighbors;
output["balance_directions"] = balance_directions;
output["uniform_mixture"] = uniform_mixture;
output["score_optimizer_maxit"] = score_optimizer_maxit;
output["coefficient_optimizer_maxit"] = coefficient_optimizer_maxit;
output["cache_score_fits"] = cache_score_fits;
output["discrete_proposal"] = discrete_proposal;
output["laplace_optimizer_maxit"] = laplace_optimizer_maxit;
output["adaptation_start"] = adapt_start_iter;
output["adaptation_end"] = adapt_end_iter;
output["numerical_engine"] = "stable_rising_factorial_v1";
output["indicator_prior_engine"] = "floating_log_gamma_v1";
output["fit_engine"] = coefficient_warm_start ?
  "deterministic_state_start_exact_context_cache_v1" : "deterministic_zero_start_exact_context_cache_v1";
output["coefficient_warm_start"] = coefficient_warm_start;
output["within_model_update"] = within_model_update;
output["within_model_attempts"] = within_model_attempts;
output["within_model_scale"] = within_model_scale;
output["within_model_fixed_probability"] = within_model_fixed_probability;
output["within_model_fixed_scale"] = within_model_fixed_scale;
output["within_fixed_proposed_by_category"] = within_fixed_proposed;
output["within_fixed_accepted_by_category"] = within_fixed_accepted;
output["within_curvature_proposed_by_category"] = within_curvature_proposed;
output["within_curvature_accepted_by_category"] = within_curvature_accepted;
output["within_model_proposed_by_category"] = within_proposed;
output["within_model_accepted_by_category"] = within_accepted;
output["alpha_adaptation_end"] = burn;
output["joint_update"]=joint_update;output["joint_scale"]=joint_scale;
output["joint_attempts"]=joint_attempts;
output["joint_proposed_by_category"]=joint_proposed;
output["joint_accepted_by_category"]=joint_accepted;
output["joint_precision_engine"]="observed_posterior_projected_dimension_scaled_v1";
output["pilot_initial_log_posterior"]=pilot_initial_log_posterior;
output["nonfinite_mh_ratios_rejected"] = dm_numeric_diagnostics::rejected_ratios.load();
output["precision_projection_calls"] = dm_numeric_diagnostics::precision_projections.load();
output["optimizer_zero_fallbacks"] = dm_numeric_diagnostics::optimizer_fallbacks.load();
output["total_seconds"] = total_seconds;
output["iterations_per_second"] = static_cast<double>(stop_iter) / total_seconds;
output["computational_engine"] = dm_runtime::use_cached ? "conditional_cached_v2" : "full_reference";
Rcpp::List continuation;
continuation["schema"]=1;continuation["completed"]=stop_iter;
continuation["config"]=checkpoint_config;continuation["elapsed_seconds"]=total_seconds;
continuation["alpha"]=alpha;
continuation["log_prior_beta"]=log_prior_beta;
continuation["log_prior_xi"]=log_prior_xi;
continuation["lambda"]=lambda;
continuation["gamma_sum_i"]=gamma_sum_i;
continuation["m_alpha"]=m_alpha;
continuation["log_g"]=log_g;
continuation["prev_log_g"]=prev_log_g;
continuation["prev_prev_log_g"]=prev_prev_log_g;
continuation["g"]=g;
continuation["prev_accept_rate_g"]=prev_accept_rate_g;
continuation["proposal_effective_window_sum"]=proposal_effective_window_sum;
continuation["proposal_effective_post_sum"]=proposal_effective_post_sum;
continuation["proposal_entropy_post_sum"]=proposal_entropy_post_sum;
continuation["proposal_max_post_sum"]=proposal_max_post_sum;
continuation["proposal_effective_this_iter"]=proposal_effective_this_iter;
continuation["proposal_entropy_this_iter"]=proposal_entropy_this_iter;
continuation["proposal_max_this_iter"]=proposal_max_this_iter;
continuation["xi"]=xi;
continuation["PP"]=PP;
continuation["within_proposed"]=within_proposed;
continuation["within_accepted"]=within_accepted;
continuation["within_fixed_proposed"]=within_fixed_proposed;
continuation["within_fixed_accepted"]=within_fixed_accepted;
continuation["within_curvature_proposed"]=within_curvature_proposed;
continuation["within_curvature_accepted"]=within_curvature_accepted;
continuation["accept_betaxi"]=accept_betaxi;
continuation["accept_move"]=accept_move;
continuation["prop_move"]=prop_move;
continuation["acc_betaxi_window"]=acc_betaxi_window;
continuation["iter_counter_g"]=iter_counter_g;
continuation["update_counter_g"]=update_counter_g;
continuation["proposal_post_count"]=proposal_post_count;
continuation["Xbeta"]=Xbeta;
continuation["alpha_mat"]=alpha_mat;
continuation["eta"]=eta;
continuation["gamma"]=gamma;
continuation["R_alpha"]=R_alpha;
continuation["Sigma_alpha"]=Sigma_alpha;
continuation["g_window_acceptance"]=g_window_acceptance;
continuation["g_window_values"]=g_window_values;
continuation["g_window_effective_neighbors"]=g_window_effective_neighbors;
continuation["log_lik"]=log_lik;
continuation["log_prior_alpha"]=log_prior_alpha;
continuation["s_alpha"]=s_alpha;
continuation["accept_alpha"]=accept_alpha;
continuation["accept_betaxi_global"]=accept_betaxi_global;
continuation["accept_move_global"]=accept_move_global;
continuation["prop_move_global"]=prop_move_global;
continuation["beta"]=beta;
continuation["rejected_ratios"]=dm_numeric_diagnostics::rejected_ratios.load();
continuation["precision_projections"]=dm_numeric_diagnostics::precision_projections.load();
continuation["optimizer_fallbacks"]=dm_numeric_diagnostics::optimizer_fallbacks.load();
output["continuation_state"]=continuation;
output["completed_iterations"]=stop_iter;
output["planned_iterations"]=runs;
output["checkpoint_engine"]="experimental_bounded_pilot_not_production_resume_v1";
return output;

}
