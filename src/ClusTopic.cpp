#include "ClusTopic.h"

using namespace Rcpp;

namespace {

void check_positive_finite(const double& x, const std::string& name)
{
  if(x <= 0.0 || !std::isfinite(x))
    throw std::runtime_error("Error in ClusTopic: " + name + " must be positive and finite");
}

void check_ClusTopic_params(const ClusTopicParams& param)
{
  check_positive_finite(param.phi, "phi");
  check_positive_finite(param.delta0, "delta0");
  check_positive_finite(param.omega, "omega");
  check_positive_finite(param.var_delta, "var_delta");

  if(param.UpdateOmega){
    check_positive_finite(param.a_omega, "a_omega");
    check_positive_finite(param.b_omega, "b_omega");
  }
}

void check_simplex_vector(const VecCol& x, const std::string& name)
{
  if(x.size() <= 0)
    throw std::runtime_error("Error in ClusTopic: " + name + " is empty");

  double tot{0.0};
  for(int i = 0; i < x.size(); i++){
    if(x(i) <= 0.0 || !std::isfinite(x(i)))
      throw std::runtime_error("Error in ClusTopic: " + name + " must have strictly positive finite entries");
    tot += x(i);
  }

  if(std::abs(tot - 1.0) > 1e-6)
    throw std::runtime_error("Error in ClusTopic: " + name + " must sum to 1");
}

void check_Lambda_star(const MatCol& Lambda_star)
{
  if(Lambda_star.rows() <= 0)
    throw std::runtime_error("Error in ClusTopic: Lambda_star has zero rows");
  if(Lambda_star.cols() <= 0)
    throw std::runtime_error("Error in ClusTopic: Lambda_star has zero cols");

  for(int k = 0; k < Lambda_star.cols(); k++){
    VecCol lambda_k = Lambda_star.col(k);
    check_simplex_vector(lambda_k, "Lambda_star column");
  }
}

MatCol sanitize_Lambda_star(const MatCol& Lambda_star)
{
  if(Lambda_star.rows() <= 0)
    throw std::runtime_error("Error in ClusTopic: Lambda_star has zero rows");
  if(Lambda_star.cols() <= 0)
    throw std::runtime_error("Error in ClusTopic: Lambda_star has zero cols");

  const double eps = 1e-12;
  MatCol out{Lambda_star};
  for(int k = 0; k < out.cols(); k++){
    double tot{0.0};
    for(int v = 0; v < out.rows(); v++){
      if(!std::isfinite(out(v,k)))
        throw std::runtime_error("Error in ClusTopic: Lambda_star contains non-finite values");
      if(out(v,k) < 0.0)
        throw std::runtime_error("Error in ClusTopic: Lambda_star contains negative values");
      if(out(v,k) == 0.0)
        out(v,k) = eps;
      tot += out(v,k);
    }
    if(tot <= 0.0)
      throw std::runtime_error("Error in ClusTopic: Lambda_star column has non-positive sum");
    out.col(k) /= tot;
  }

  return out;
}

void check_zeta_vector(const VecCol& zeta, const unsigned int& V)
{
  if(zeta.size() != V)
    throw std::runtime_error("Error in ClusTopic: zeta has incompatible size");

  for(int v = 0; v < zeta.size(); v++){
    if(zeta(v) <= 0.0 || !std::isfinite(zeta(v)))
      throw std::runtime_error("Error in ClusTopic: zeta must have strictly positive finite entries");
  }
}

double log_sum_exp(const VecCol& log_w)
{
  if(log_w.size() <= 0)
    throw std::runtime_error("Error in ClusTopic: empty log weights");

  for(int i = 0; i < log_w.size(); i++){
    if(std::isnan(log_w(i)))
      throw std::runtime_error("Error in ClusTopic: log weights contain NaN values");
  }

  double max_log_w = log_w.maxCoeff();
  if(!std::isfinite(max_log_w))
    throw std::runtime_error("Error in ClusTopic: all log weights are non-finite");

  double tot{0.0};
  for(int i = 0; i < log_w.size(); i++)
    tot += std::exp(log_w(i) - max_log_w);

  return max_log_w + std::log(tot);
}

VecCol normalize_log_weights(const VecCol& log_w)
{
  VecCol prob{VecCol::Zero(log_w.size())};
  double lse = log_sum_exp(log_w);

  for(int i = 0; i < log_w.size(); i++)
    prob(i) = std::exp(log_w(i) - lse);

  return prob;
}

VecCol sample_dirichlet_symmetric(sample::GSL_RNG const & engine, const unsigned int& V, const double& delta0)
{
  sample::rgamma rgamma;
  VecCol out{VecCol::Zero(V)};

  double tot{0.0};
  for(unsigned int v = 0; v < V; v++){
    out(v) = rgamma(engine, delta0, 1.0);
    if(out(v) <= 0.0 || !std::isfinite(out(v)))
      out(v) = std::numeric_limits<double>::min();
    tot += out(v);
  }

  if(tot <= 0.0 || !std::isfinite(tot)){
    out = VecCol::Constant(V, 1.0/(double)V);
    return out;
  }

  out /= tot;
  return out;
}

VecCol simplex_to_logratio(const VecCol& delta)
{
  // Map a simplex vector delta = (delta_1,...,delta_V) to V-1 unconstrained
  // log-ratios, using the last component as reference:
  // eta_v = log(delta_v / delta_V), v = 1,...,V-1.
  // This parametrization lets the MH proposal move delta in Euclidean space.
  if(delta.size() == 1)
    return VecCol::Zero(0);

  VecCol eta{VecCol::Zero(delta.size() - 1)};
  const double ref = delta(delta.size() - 1);
  for(int v = 0; v < eta.size(); v++)
    eta(v) = std::log(delta(v)) - std::log(ref);

  return eta;
}

VecCol logratio_to_simplex(const VecCol& eta)
{
  // Inverse map of simplex_to_logratio. Given eta_v = log(delta_v/delta_V),
  // recover the simplex vector:
  // delta_v = exp(eta_v)/(1 + sum_j exp(eta_j)), v = 1,...,V-1,
  // delta_V = 1/(1 + sum_j exp(eta_j)).
  // The max_eta shift is only for numerical stability.
  const int V = eta.size() + 1;
  if(V == 1)
    return VecCol::Constant(1, 1.0);

  double max_eta = 0.0;
  for(int v = 0; v < eta.size(); v++)
    max_eta = std::max(max_eta, eta(v));

  double denom = std::exp(-max_eta);
  for(int v = 0; v < eta.size(); v++)
    denom += std::exp(eta(v) - max_eta);

  VecCol delta{VecCol::Zero(V)};
  for(int v = 0; v < eta.size(); v++)
    delta(v) = std::exp(eta(v) - max_eta)/denom;
  delta(V - 1) = std::exp(-max_eta)/denom;

  return delta;
}

double log_jacobian_delta(const VecCol& delta)
{
  double res{0.0};
  for(int v = 0; v < delta.size(); v++)
    res += std::log(delta(v));

  return res;
}

} // anonymous namespace


double log_ClusTopic_dirichlet_density(const VecCol& lambda, const VecCol& zeta)
{
  if(lambda.size() != zeta.size())
    throw std::runtime_error("Error in log_ClusTopic_dirichlet_density: lambda and zeta have incompatible sizes");

  check_simplex_vector(lambda, "lambda");
  check_zeta_vector(zeta, lambda.size());

  double alpha0 = zeta.sum();
  double res = std::lgamma(alpha0);
  for(int v = 0; v < lambda.size(); v++){
    res -= std::lgamma(zeta(v));
    res += (zeta(v) - 1.0)*std::log(lambda(v));
  }

  return res;
}


double log_ClusTopic_allocated_zeta_full_conditional(const VecCol& A_m, const unsigned int& n_m,
                                                     const double& phi, const VecCol& delta,
                                                     const double& delta0)
{
  if(n_m == 0)
    throw std::runtime_error("Error in log_ClusTopic_allocated_zeta_full_conditional: n_m must be positive");
  if(A_m.size() != delta.size())
    throw std::runtime_error("Error in log_ClusTopic_allocated_zeta_full_conditional: A_m and delta have incompatible sizes");

  check_positive_finite(phi, "phi");
  check_positive_finite(delta0, "delta0");
  check_simplex_vector(delta, "delta");

  double res = (delta0 - 1.0)*delta.array().log().sum();
  res += (double)n_m*std::lgamma(phi);

  for(int v = 0; v < delta.size(); v++)
    res -= (double)n_m*std::lgamma(phi*delta(v));

  res += phi*delta.dot(A_m);
  return res;
}


VecUnsCol compute_ClusTopic_cluster_sizes(const VecIntCol& c, const unsigned int& M)
{
  if(M <= 0)
    throw std::runtime_error("Error in compute_ClusTopic_cluster_sizes: M must be positive");

  VecUnsCol cluster_size{VecUnsCol::Zero(M)};
  for(int k = 0; k < c.size(); k++){
    if(c(k) < 0 || c(k) >= (int)M)
      throw std::runtime_error("Error in compute_ClusTopic_cluster_sizes: c contains invalid labels");
    cluster_size(c(k))++;
  }

  return cluster_size;
}


MatCol compute_ClusTopic_A(const MatCol& Lambda_star, const VecIntCol& c, const unsigned int& M)
{
  check_Lambda_star(Lambda_star);
  if(c.size() != Lambda_star.cols())
    throw std::runtime_error("Error in compute_ClusTopic_A: c has incompatible size");

  const int V = Lambda_star.rows();
  const int K = Lambda_star.cols();
  MatCol A{MatCol::Zero(V, M)};

  for(int k = 0; k < K; k++){
    if(c(k) < 0 || c(k) >= (int)M)
      throw std::runtime_error("Error in compute_ClusTopic_A: c contains invalid labels");
    A.col(c(k)) += Lambda_star.col(k).array().log().matrix();
  }

  return A;
}


ClusTopicZetaDraw split_ClusTopic_zeta(const VecCol& zeta, const double& phi)
{
  if(zeta.size() <= 0)
    throw std::runtime_error("Error in split_ClusTopic_zeta: zeta is empty");

  check_zeta_vector(zeta, zeta.size());
  check_positive_finite(phi, "phi");

  ClusTopicZetaDraw out;
  out.phi = phi;
  out.delta = zeta/zeta.sum();
  check_simplex_vector(out.delta, "delta");
  out.zeta = out.phi*out.delta;
  return out;
}


ClusTopicZetaDraw sample_ClusTopic_prior_zeta(sample::GSL_RNG const & engine, const unsigned int& V,
                                              const double& phi, const double& delta0)
{
  if(V <= 0)
    throw std::runtime_error("Error in sample_ClusTopic_prior_zeta: V must be positive");
  check_positive_finite(phi, "phi");
  check_positive_finite(delta0, "delta0");

  ClusTopicZetaDraw out;
  out.phi = phi;
  out.delta = sample_dirichlet_symmetric(engine, V, delta0);
  out.zeta = out.phi*out.delta;
  return out;
}


unsigned int sample_ClusTopic_Mstar(sample::GSL_RNG const & engine, const unsigned int& M,
                                    const unsigned int& K, const double& omega,
                                    const unsigned int& mstar_max,
                                    VecCol& mstar_prob, VecCol& log_mstar_prob)
{
  if(M <= 0)
    throw std::runtime_error("Error in sample_ClusTopic_Mstar: M must be positive");
  if(K <= 0)
    throw std::runtime_error("Error in sample_ClusTopic_Mstar: K must be positive");
  check_positive_finite(omega, "omega");

  const unsigned int support_size = mstar_max + 1;
  log_mstar_prob = VecCol::Zero(support_size);

  for(unsigned int mstar = 0; mstar < support_size; mstar++){
    double lp = (double)mstar*std::log(omega) - std::lgamma((double)mstar + 1.0);
    if(K > 1)
      lp -= ((double)K - 1.0)*std::log((double)M + (double)mstar);
    log_mstar_prob(mstar) = lp;
  }

  mstar_prob = normalize_log_weights(log_mstar_prob);

  sample::sample_index rsample;
  return rsample(engine, mstar_prob);
}


ClusTopicZetaMH sample_ClusTopic_allocated_zeta(sample::GSL_RNG const & engine, const VecCol& zeta_old,
                                                const VecCol& A_m, const unsigned int& n_m,
                                                const double& phi, const double& delta0,
                                                const double& var_delta)
{
  check_positive_finite(phi, "phi");
  check_positive_finite(var_delta, "var_delta");

  sample::rnorm rnorm;
  sample::runif runif;

  ClusTopicZetaDraw old_draw = split_ClusTopic_zeta(zeta_old, phi);
  VecCol eta_old = simplex_to_logratio(old_draw.delta);

  VecCol eta_prime{eta_old};
  for(int v = 0; v < eta_prime.size(); v++)
    eta_prime(v) = rnorm(engine, eta_old(v), std::sqrt(var_delta));

  VecCol delta_prime = logratio_to_simplex(eta_prime);

  double log_old = log_ClusTopic_allocated_zeta_full_conditional(A_m, n_m, phi, old_draw.delta,
                                                                 delta0);
  double log_prime = log_ClusTopic_allocated_zeta_full_conditional(A_m, n_m, phi, delta_prime,
                                                                   delta0);

  log_old += log_jacobian_delta(old_draw.delta);
  log_prime += log_jacobian_delta(delta_prime);

  ClusTopicZetaMH out;
  out.log_acc = log_prime - log_old;

  double log_u = std::log(runif(engine));
  if(log_u < std::min(0.0, out.log_acc)){
    out.accepted = true;
    out.phi = phi;
    out.delta = delta_prime;
    out.zeta = phi*delta_prime;
  }
  else{
    out.accepted = false;
    out.phi = phi;
    out.delta = old_draw.delta;
    out.zeta = old_draw.zeta;
  }

  return out;
}


ClusTopicUpdate sample_ClusTopic_partition(sample::GSL_RNG const & engine, const MatCol& Lambda_star,
                                           const std::vector<VecCol>& Zeta_old,
                                           const ClusTopicParams& param)
{
  // 1) Basic validation. Lambda_star is V x K, with each column a simplex
  // vector; Zeta_old contains all currently available component centers
  // (allocated and empty) before this clustering update.
  check_ClusTopic_params(param);
  MatCol Lambda_star_work = sanitize_Lambda_star(Lambda_star);
  check_Lambda_star(Lambda_star_work);

  const unsigned int V = Lambda_star_work.rows();
  const unsigned int K = Lambda_star_work.cols();
  const unsigned int T_old = Zeta_old.size();
  if(T_old <= 0)
    throw std::runtime_error("Error in sample_ClusTopic_partition: Zeta_old is empty");

  for(unsigned int m = 0; m < T_old; m++)
    check_zeta_vector(Zeta_old[m], V);

  // 2) Allocation update. For each topic vector lambda_k, compute the
  // Dirichlet log-density under every old center zeta_m, normalize the log
  // weights, and sample the raw allocation label.
  sample::sample_index rsample;
  ClusTopicUpdate out;
  out.aux.c_raw = VecIntCol::Zero(K);
  out.aux.allocation_prob = MatCol::Zero(K, T_old);

  for(unsigned int k = 0; k < K; k++){
    VecCol log_prob{VecCol::Zero(T_old)};
    VecCol lambda_k = Lambda_star_work.col(k);

    for(unsigned int m = 0; m < T_old; m++)
      log_prob(m) = log_ClusTopic_dirichlet_density(lambda_k, Zeta_old[m]);

    VecCol prob = normalize_log_weights(log_prob);
    out.aux.allocation_prob.row(k) = prob.transpose();
    out.aux.c_raw(k) = (int)rsample(engine, prob);
  }

  // 3) Remove empty labels created by the allocation step. The resulting
  // partition uses compact labels 0,...,M-1, and Zeta_alloc stores the old
  // centers associated with these allocated clusters.
  std::vector<int> relabel(T_old, -1);
  std::vector<VecCol> Zeta_alloc;
  out.c = VecIntCol::Zero(K);

  for(unsigned int k = 0; k < K; k++){
    const int old_label = out.aux.c_raw(k);
    if(relabel[old_label] < 0){
      relabel[old_label] = (int)Zeta_alloc.size();
      Zeta_alloc.push_back(Zeta_old[old_label]);
    }
    out.c(k) = relabel[old_label];
  }

  out.M = Zeta_alloc.size();

  // 4) Build sufficient statistics for the allocated centers. cluster_size[m]
  // is n_m, while A.col(m) stores A_{m,v} = sum_{k:c_k=m} log(lambda_{v,k}).
  out.aux.cluster_size = compute_ClusTopic_cluster_sizes(out.c, out.M);
  out.aux.A = compute_ClusTopic_A(Lambda_star_work, out.c, out.M);

  // 5) Update the number of empty components. The infinite support in the
  // algorithm is truncated to 0,...,param.mstar_max.
  out.Mstar = sample_ClusTopic_Mstar(engine, out.M, K, param.omega, param.mstar_max,
                                     out.aux.mstar_prob, out.aux.log_mstar_prob);

  // 6) Prepare output containers. Centers are returned with allocated centers
  // first, followed by Mstar empty centers drawn from the prior.
  const unsigned int T_new = out.M + out.Mstar;
  out.Zeta.reserve(T_new);
  out.aux.delta.reserve(T_new);
  out.aux.phi = param.phi;
  out.aux.log_acc_zeta = VecCol::Zero(out.M);
  out.aux.accept_zeta = VecUnsCol::Zero(out.M);

  // 7) Update allocated centers. If requested, each allocated zeta_m is moved
  // with an MH step on log-ratios(delta_m), keeping the common phi fixed.
  // Otherwise it is simply decomposed into delta_m and rescaled by phi.
  for(unsigned int m = 0; m < out.M; m++){
    if(param.UpdateZeta){
      ClusTopicZetaMH mh = sample_ClusTopic_allocated_zeta(engine, Zeta_alloc[m], out.aux.A.col(m),
                                                           out.aux.cluster_size(m),
                                                           param.phi, param.delta0,
                                                           param.var_delta);
      out.Zeta.push_back(mh.zeta);
      out.aux.delta.push_back(mh.delta);
      out.aux.log_acc_zeta(m) = mh.log_acc;
      out.aux.accept_zeta(m) = mh.accepted ? 1 : 0;
    }
    else{
      ClusTopicZetaDraw draw = split_ClusTopic_zeta(Zeta_alloc[m], param.phi);
      out.Zeta.push_back(draw.zeta);
      out.aux.delta.push_back(draw.delta);
    }
  }

  // 8) Refresh empty component parameters from the prior, since empty clusters
  // have no likelihood contribution.
  for(unsigned int mstar = 0; mstar < out.Mstar; mstar++){
    ClusTopicZetaDraw draw = sample_ClusTopic_prior_zeta(engine, V, param.phi, param.delta0);
    out.Zeta.push_back(draw.zeta);
    out.aux.delta.push_back(draw.delta);
  }

  // 9) Optional hyperprior update for omega. The shape/rate are stored even
  // when omega is not sampled, so the caller can inspect the full conditional.
  out.aux.omega_shape = param.a_omega + (double)out.M + (double)out.Mstar - 1.0;
  out.aux.omega_rate = param.b_omega + 1.0;
  out.omega = param.omega;
  if(param.UpdateOmega){
    sample::rgamma rgamma;
    out.omega = rgamma(engine, out.aux.omega_shape, 1.0/out.aux.omega_rate);
  }

  return out;
}


ClusTopicUpdate sample_ClusTopic_partition_fixedM(sample::GSL_RNG const & engine, const MatCol& Lambda_star,
                                                  const std::vector<VecCol>& Zeta_old,
                                                  const ClusTopicParams& param)
{
  // Fixed-M version used by the DTM random-center sampler. The number of
  // centers is kept equal to Zeta_old.size(); empty clusters are allowed, but
  // Mstar is bypassed and set to zero.
  check_ClusTopic_params(param);
  MatCol Lambda_star_work = sanitize_Lambda_star(Lambda_star);
  check_Lambda_star(Lambda_star_work);

  const unsigned int V = Lambda_star_work.rows();
  const unsigned int K = Lambda_star_work.cols();
  const unsigned int M = Zeta_old.size();
  if(M <= 0)
    throw std::runtime_error("Error in sample_ClusTopic_partition_fixedM: Zeta_old is empty");

  for(unsigned int m = 0; m < M; m++)
    check_zeta_vector(Zeta_old[m], V);

  sample::sample_index rsample;
  ClusTopicUpdate out;
  out.M = M;
  out.Mstar = 0;
  out.omega = param.omega;
  out.c = VecIntCol::Zero(K);
  out.aux.c_raw = VecIntCol::Zero(K);
  out.aux.allocation_prob = MatCol::Zero(K, M);

  // Allocate each Lambda_star column to one of the existing centers. Labels
  // are not relabelled/compacted, so a center can have cluster_size[m] = 0.
  for(unsigned int k = 0; k < K; k++){
    VecCol log_prob{VecCol::Zero(M)};
    VecCol lambda_k = Lambda_star_work.col(k);

    for(unsigned int m = 0; m < M; m++)
      log_prob(m) = log_ClusTopic_dirichlet_density(lambda_k, Zeta_old[m]);

    VecCol prob = normalize_log_weights(log_prob);
    const int label = (int)rsample(engine, prob);
    out.aux.allocation_prob.row(k) = prob.transpose();
    out.aux.c_raw(k) = label;
    out.c(k) = label;
  }

  out.aux.cluster_size = compute_ClusTopic_cluster_sizes(out.c, out.M);
  out.aux.A = compute_ClusTopic_A(Lambda_star_work, out.c, out.M);

  // Mstar is intentionally not sampled in this fixed-M regime. We still fill
  // the auxiliary probability vectors with the degenerate distribution at 0.
  out.aux.mstar_prob = VecCol::Constant(1, 1.0);
  out.aux.log_mstar_prob = VecCol::Constant(1, 0.0);

  out.Zeta.reserve(M);
  out.aux.delta.reserve(M);
  out.aux.phi = param.phi;
  out.aux.log_acc_zeta = VecCol::Zero(M);
  out.aux.accept_zeta = VecUnsCol::Zero(M);

  for(unsigned int m = 0; m < M; m++){
    if(param.UpdateZeta && out.aux.cluster_size(m) > 0){
      ClusTopicZetaMH mh = sample_ClusTopic_allocated_zeta(engine, Zeta_old[m], out.aux.A.col(m),
                                                           out.aux.cluster_size(m),
                                                           param.phi, param.delta0,
                                                           param.var_delta);
      out.Zeta.push_back(mh.zeta);
      out.aux.delta.push_back(mh.delta);
      out.aux.log_acc_zeta(m) = mh.log_acc;
      out.aux.accept_zeta(m) = mh.accepted ? 1 : 0;
    }
    else if(param.UpdateZeta){
      // Empty fixed centers have no likelihood contribution.
      ClusTopicZetaDraw draw = sample_ClusTopic_prior_zeta(engine, V, param.phi, param.delta0);
      out.Zeta.push_back(draw.zeta);
      out.aux.delta.push_back(draw.delta);
    }
    else{
      ClusTopicZetaDraw draw = split_ClusTopic_zeta(Zeta_old[m], param.phi);
      out.Zeta.push_back(draw.zeta);
      out.aux.delta.push_back(draw.delta);
    }
  }

  out.aux.omega_shape = param.a_omega + (double)M - 1.0;
  out.aux.omega_rate = param.b_omega + 1.0;
  if(param.UpdateOmega){
    sample::rgamma rgamma;
    out.omega = rgamma(engine, out.aux.omega_shape, 1.0/out.aux.omega_rate);
  }

  return out;
}
