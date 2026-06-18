#include "DTMcenters.h"

using namespace Rcpp;

namespace {

const double PROCESS_PHI = 1.0;

void check_positive_finite_centers(const double& x, const std::string& name)
{
  if(x <= 0.0 || !std::isfinite(x))
    throw std::runtime_error("Error in DTM centers: " + name + " must be positive and finite");
}

void check_zeta_for_lambda_update(const VecCol& zeta, const unsigned int& V)
{
  if(zeta.size() != V)
    throw std::runtime_error("Error in DTM centers: zeta has incompatible size");

  for(int v = 0; v < zeta.size(); v++){
    if(zeta(v) <= 0.0 || !std::isfinite(zeta(v)))
      throw std::runtime_error("Error in DTM centers: zeta must contain positive finite values");
  }
}

void check_center_state_shapes(const DTMcenterState& state, const unsigned int& H,
                               const unsigned int& V, const unsigned int& Ttot)
{
  if(state.M == 0)
    throw std::runtime_error("Error in DTM centers: M must be positive");
  if(state.Xi.size() != state.M || state.S.size() != state.M || state.U.size() != state.M ||
     state.Lambda.size() != state.M || state.Dl.size() != state.M || state.N.size() != state.M)
    throw std::runtime_error("Error in DTM centers: process-state lists must all have length M");
  if(state.Zeta.size() < state.M)
    throw std::runtime_error("Error in DTM centers: Zeta must contain at least M centers");

  for(unsigned int m = 0; m < state.M; m++){
    if(state.Xi[m].rows() != H || state.Xi[m].cols() != Ttot)
      throw std::runtime_error("Error in DTM centers: invalid Xi dimensions");
    if(state.S[m].rows() != H || state.S[m].cols() != Ttot)
      throw std::runtime_error("Error in DTM centers: invalid S dimensions");
    if(state.U[m].rows() != H || state.U[m].cols() != Ttot)
      throw std::runtime_error("Error in DTM centers: invalid U dimensions");
    if(state.N[m].rows() != Ttot || state.N[m].cols() != H)
      throw std::runtime_error("Error in DTM centers: invalid N dimensions");
    if(state.Lambda[m].size() != H || state.Dl[m].size() != H)
      throw std::runtime_error("Error in DTM centers: invalid Lambda/Dl list length");

    check_zeta_for_lambda_update(state.Zeta[m], V);
    for(unsigned int l = 0; l < H; l++){
      if(state.Lambda[m][l].rows() != V || state.Lambda[m][l].cols() != Ttot)
        throw std::runtime_error("Error in DTM centers: invalid Lambda dimensions");
      if(state.Dl[m][l].rows() != V || state.Dl[m][l].cols() != Ttot)
        throw std::runtime_error("Error in DTM centers: invalid Dl dimensions");
    }
  }
}

std::vector<VecCol> allocated_zeta(const std::vector<VecCol>& Zeta, const unsigned int& M)
{
  if(Zeta.size() < M)
    throw std::runtime_error("Error in allocated_zeta: Zeta.size() < M");

  std::vector<VecCol> out(M);
  for(unsigned int m = 0; m < M; m++)
    out[m] = Zeta[m];
  return out;
}

std::vector<MatCol> initialize_Lambda_center(sample::GSL_RNG const & engine,
                                             const unsigned int& H,
                                             const unsigned int& V,
                                             const unsigned int& Ttot,
                                             const VecCol& zeta)
{
  MatIntCol Xi_empty{MatIntCol::Zero(H, Ttot)};
  Xi_empty.row(0).setOnes();

  std::vector<MatUnsCol> Dl_empty(H, MatUnsCol::Zero(V, Ttot));
  return sample_Lambda_itlm(engine, Dl_empty, Xi_empty, zeta);
}

MatIntCol initialize_Xi_center(const unsigned int& H, const unsigned int& Ttot)
{
  MatIntCol Xi_new{MatIntCol::Zero(H, Ttot)};
  Xi_new.row(0).setOnes();
  return Xi_new;
}

} // anonymous namespace


double compute_t_sigma_gamma_centers(const double& gamma, const double& sigma,
                                     const unsigned int& H, const unsigned int& M)
{
  check_positive_finite_centers(gamma, "gamma");
  if(sigma <= 0.0 || sigma >= 1.0 || !std::isfinite(sigma))
    throw std::runtime_error("Error in compute_t_sigma_gamma_centers: sigma must be in (0,1)");
  if(H == 0 || M == 0)
    throw std::runtime_error("Error in compute_t_sigma_gamma_centers: H and M must be positive");

  const double gamma_scaled = gamma/(double)M;
  return std::exp((std::log(sigma*(double)H) - std::log(gamma_scaled))/sigma);
}


MatIntCol sample_Xi_tl_centers(sample::GSL_RNG const & engine,
                               const MatIntCol& Xi_old,
                               const MatCol& S,
                               const MatUnsCol& N_tl,
                               const double& phi,
                               const double& sigma,
                               const double& b,
                               const double& t_sigma_gamma)
{
  // Same MH update as sample_Xi_tl(), but suitable for the fixed-center
  // sampler: a single center/process is allowed to be inactive at a time t.
  // Therefore we do not enforce positive column sums for this individual Xi_m.
  sample::rpoisson rpoisson;
  sample::runif runif;

  const int H = S.rows();
  const int Ttot = S.cols();
  if(H <= 0)
    throw std::runtime_error("Error in sample_Xi_tl_centers: H must be positive");
  if(Xi_old.cols() != Ttot || Xi_old.rows() != H)
    throw std::runtime_error("Error in sample_Xi_tl_centers: Xi_old has incompatible dimensions");
  if(N_tl.rows() != Ttot || N_tl.cols() != H)
    throw std::runtime_error("Error in sample_Xi_tl_centers: N_tl has incompatible dimensions");
  if(phi <= 0.0 || b <= 0.0)
    throw std::runtime_error("Error in sample_Xi_tl_centers: phi and beta must be positive");
  if(sigma < 0.0 || sigma >= 1.0)
    throw std::runtime_error("Error in sample_Xi_tl_centers: sigma is out of range");
  if(t_sigma_gamma <= 0.0 || !std::isfinite(t_sigma_gamma))
    throw std::runtime_error("Error in sample_Xi_tl_centers: invalid t_sigma_gamma");

  MatIntCol Xi{MatIntCol::Zero(H,Ttot)};
  const double b_phi = b + phi;
  const double b_phi_t = b + phi + t_sigma_gamma;
  const double diff_log = std::log(b_phi) - std::log(b_phi_t);

  auto logNormConst = [b_phi,b_phi_t,sigma,diff_log](int c){
    double res{0.0};
    if(c == 0){
      res += std::lgamma(1.0 - sigma) - std::log(sigma) +
             std::log(gsl_expm1(-sigma*diff_log));
    }
    else if(c >= 1){
      res += std::lgamma((double)c - sigma) +
             std::log(-gsl_expm1(((double)c - sigma)*diff_log));
    }
    else{
      throw std::runtime_error("Error in sample_Xi_tl_centers: c must be non-negative");
    }
    return res;
  };

  for(int l = 0; l < H; l++){
    for(int t = 0; t < Ttot; t++){
      if(S(l,t) <= 0.0)
        throw std::runtime_error("Error in sample_Xi_tl_centers: S_tl must be positive");

      const int xi_prime = rpoisson(engine, phi*S(l,t));
      const int xi_tl = Xi_old(l,t);
      double pacc{0.0};

      if(!(N_tl(t,l) > 0 && xi_prime == 0)){
        const int diff_xi = xi_prime - xi_tl;
        double log_R_tl = -(double)diff_xi;

        if(xi_tl > 0 && xi_prime > 0){
          log_R_tl += N_tl(t,l) *
            (std::log((double)xi_prime) - std::log((double)xi_tl));
        }

        if(t < Ttot - 1){
          log_R_tl += (double)diff_xi *
            (std::log(S(l,t+1)) + std::log(b_phi) - b_phi);
          log_R_tl += logNormConst(xi_tl) - logNormConst(xi_prime);
        }

        if(!std::isfinite(log_R_tl))
          throw std::runtime_error("Error in sample_Xi_tl_centers: invalid acceptance probability");

        pacc = std::exp(std::min(0.0, log_R_tl));
      }

      if(runif(engine) < pacc)
        Xi(l,t) = xi_prime;
      else
        Xi(l,t) = xi_tl;
    }
  }

  return Xi;
}


std::vector<MatCol> sample_Lambda_itlm(sample::GSL_RNG const & engine,
                                       const std::vector<MatUnsCol>& D_itl,
                                       const MatIntCol& Xi,
                                       const VecCol& zeta)
{
  sample::rgamma rgamma;
  sample::sample_index rsample;

  const int H = Xi.rows();
  const int Ttot = Xi.cols();
  if(H <= 0)
    throw std::runtime_error("Error in sample_Lambda_itlm: H must be positive");
  if(D_itl.size() != H)
    throw std::runtime_error("Error in sample_Lambda_itlm: D_itl has incompatible length");
  if(D_itl[0].cols() != Ttot)
    throw std::runtime_error("Error in sample_Lambda_itlm: D_itl has incompatible time dimension");

  const int V = D_itl[0].rows();
  check_zeta_for_lambda_update(zeta, V);

  std::vector<MatCol> Lambda_itl(H, MatCol::Zero(V, Ttot));

  for(int l = 0; l < H; l++){
    std::vector<int> idx_born;
    std::vector<int> idx_surv;
    std::vector<int> idx_noact;

    VecIntCol Xi_l = Xi.row(l);
    find_indices(Xi.row(l), idx_born, idx_surv, idx_noact);

    for(int t : idx_born){
      VecCol temp{VecCol::Zero(V)};
      for(int i = 0; i < V; i++){
        double shape = zeta(i);
        int s = t;
        bool flag = true;
        while(flag && s < Ttot){
          shape += (double)D_itl[l](i, s);
          if(s + 1 >= Ttot || Xi_l[s + 1] == 0)
            flag = false;
          s++;
        }
        temp(i) = rgamma(engine, shape, 1.0);
        if(temp(i) <= 0.0 || !std::isfinite(temp(i)))
          temp(i) = std::numeric_limits<double>::min();
      }
      if(temp.sum() <= 0.0)
        temp(rsample(engine, temp.size())) = 1.0;
      temp /= temp.sum();
      Lambda_itl[l].col(t) = temp;
    }

    for(int t : idx_surv){
      if(t == 0)
        throw std::runtime_error("Error in sample_Lambda_itlm: survivor at first time");
      Lambda_itl[l].col(t) = Lambda_itl[l].col(t - 1);
    }

    for(int t : idx_noact){
      VecCol temp{VecCol::Zero(V)};
      for(int i = 0; i < V; i++){
        temp(i) = rgamma(engine, zeta(i), 1.0);
        if(temp(i) <= 0.0 || !std::isfinite(temp(i)))
          temp(i) = std::numeric_limits<double>::min();
      }
      if(temp.sum() <= 0.0)
        temp(rsample(engine, temp.size())) = 1.0;
      temp /= temp.sum();
      Lambda_itl[l].col(t) = temp;
    }
  }

  return Lambda_itl;
}


std::pair<std::vector<std::vector<MatUnsCol>>, std::vector<MatUnsCol>>
sample_Ditlm(sample::GSL_RNG const & engine,
             const std::vector<std::vector<MatCol>>& Lambda,
             const std::vector<MatIntCol>& Xi,
             const MatIntCol& D)
{
  sample::rmultinomial<VecUnsCol> rmultinomial;

  const unsigned int M = Xi.size();
  if(M == 0)
    throw std::runtime_error("Error in sample_Ditlm: M must be positive");
  if(Lambda.size() != M)
    throw std::runtime_error("Error in sample_Ditlm: Lambda.size() != M");

  const unsigned int H = Xi[0].rows();
  const unsigned int Ttot = Xi[0].cols();
  const unsigned int V = D.rows();
  if(D.cols() != Ttot)
    throw std::runtime_error("Error in sample_Ditlm: D has incompatible time dimension");

  std::vector<std::vector<MatUnsCol>> Dl(M, std::vector<MatUnsCol>(H, MatUnsCol::Zero(V, Ttot)));
  std::vector<MatUnsCol> N(M, MatUnsCol::Zero(Ttot, H));

  for(unsigned int m = 0; m < M; m++){
    if(Xi[m].rows() != H || Xi[m].cols() != Ttot)
      throw std::runtime_error("Error in sample_Ditlm: Xi dimensions are not consistent");
    if(Lambda[m].size() != H)
      throw std::runtime_error("Error in sample_Ditlm: Lambda[m].size() != H");
  }

  for(unsigned int i = 0; i < V; i++){
    for(unsigned int t = 0; t < Ttot; t++){
      VecCol weights{VecCol::Zero(M*H)};
      for(unsigned int m = 0; m < M; m++){
        for(unsigned int l = 0; l < H; l++){
          if(Lambda[m][l].rows() != V || Lambda[m][l].cols() != Ttot)
            throw std::runtime_error("Error in sample_Ditlm: invalid Lambda dimensions");
          weights(m*H + l) = Lambda[m][l](i, t) * (double)Xi[m](l, t);
        }
      }

      const double sum_w = weights.sum();
      if(sum_w <= 0.0){
        if(D(i, t) > 0)
          throw std::runtime_error("Error in sample_Ditlm: positive count with zero total weight");
        continue;
      }

      weights /= sum_w;
      VecUnsCol temp = rmultinomial(engine, D(i, t), weights);
      for(unsigned int m = 0; m < M; m++){
        for(unsigned int l = 0; l < H; l++)
          Dl[m][l](i, t) = temp(m*H + l);
      }
    }
  }

  for(unsigned int m = 0; m < M; m++){
    for(unsigned int l = 0; l < H; l++)
      N[m].col(l) = Dl[m][l].colwise().sum();
  }

  return std::make_pair(Dl, N);
}


DTMcenterLambdaStar build_Lambda_star_centers(const std::vector<std::vector<MatCol>>& Lambda,
                                              const std::vector<MatIntCol>& Xi)
{
  const unsigned int M = Xi.size();
  if(M == 0)
    throw std::runtime_error("Error in build_Lambda_star_centers: M must be positive");

  const unsigned int H = Xi[0].rows();
  const unsigned int V = Lambda[0][0].rows();
  std::vector<std::vector<VecCol>> tmp(M);
  std::vector<int> center_id;
  std::vector<int> row_id;
  std::vector<int> born_time;

  for(unsigned int m = 0; m < M; m++){
    if(Lambda[m].size() != H)
      throw std::runtime_error("Error in build_Lambda_star_centers: invalid Lambda length");

    for(unsigned int l = 0; l < H; l++){
      std::vector<int> idx_born;
      std::vector<int> idx_surv;
      std::vector<int> idx_noact;
      find_indices(Xi[m].row(l), idx_born, idx_surv, idx_noact);

      for(int t : idx_born){
        tmp[m].push_back(Lambda[m][l].col(t));
        center_id.push_back((int)m);
        row_id.push_back((int)l);
        born_time.push_back(t);
      }
    }
  }

  unsigned int K{0};
  for(unsigned int m = 0; m < M; m++)
    K += tmp[m].size();

  DTMcenterLambdaStar out;
  out.by_center.resize(M);
  out.all = MatCol::Zero(V, K);
  out.center_id = VecIntCol::Zero(K);
  out.row_id = VecIntCol::Zero(K);
  out.born_time = VecIntCol::Zero(K);

  unsigned int pos{0};
  for(unsigned int m = 0; m < M; m++){
    out.by_center[m] = MatCol::Zero(V, tmp[m].size());
    for(unsigned int k = 0; k < tmp[m].size(); k++){
      out.by_center[m].col(k) = tmp[m][k];
      out.all.col(pos) = tmp[m][k];
      out.center_id(pos) = center_id[pos];
      out.row_id(pos) = row_id[pos];
      out.born_time(pos) = born_time[pos];
      pos++;
    }
  }

  return out;
}


void resize_DTMcenter_state_after_clustering(sample::GSL_RNG const & engine,
                                             DTMcenterState& state,
                                             const unsigned int& M_new,
                                             const unsigned int& H,
                                             const unsigned int& V,
                                             const unsigned int& Ttot,
                                             const double& gamma,
                                             const double& sigma,
                                             const double& beta)
{
  if(M_new == 0)
    throw std::runtime_error("Error in resize_DTMcenter_state_after_clustering: M_new must be positive");
  if(state.Zeta.size() < M_new)
    throw std::runtime_error("Error in resize_DTMcenter_state_after_clustering: not enough centers");

  const unsigned int M_old = state.M;
  if(M_new < M_old){
    state.Xi.resize(M_new);
    state.S.resize(M_new);
    state.U.resize(M_new);
    state.Lambda.resize(M_new);
    state.Dl.resize(M_new);
    state.N.resize(M_new);
  }
  else if(M_new > M_old){
    const double t_sigma_gamma = compute_t_sigma_gamma_centers(gamma, sigma, H, M_new);
    for(unsigned int m = M_old; m < M_new; m++){
      MatIntCol Xi_new = initialize_Xi_center(H, Ttot);
      MatCol U_new = MatCol::Zero(H, Ttot);
      MatCol S_new = sample_Stl(engine, Xi_new, U_new, PROCESS_PHI, sigma, beta);
      U_new = sample_Utl(engine, S_new, t_sigma_gamma);
      std::vector<MatCol> Lambda_new = initialize_Lambda_center(engine, H, V, Ttot, state.Zeta[m]);
      std::vector<MatUnsCol> Dl_new(H, MatUnsCol::Zero(V, Ttot));
      MatUnsCol N_new = MatUnsCol::Zero(Ttot, H);

      state.Xi.push_back(Xi_new);
      state.S.push_back(S_new);
      state.U.push_back(U_new);
      state.Lambda.push_back(Lambda_new);
      state.Dl.push_back(Dl_new);
      state.N.push_back(N_new);
    }
  }

  state.M = M_new;
}


Rcpp::List GibbsSampler_DTM_centers_c_core(const int& niter, const int& nburn, const int& thin,
                                           const MatIntCol& D, const int& H, const int& V,
                                           const int& Ttot, const DTMcenterParams& param,
                                           DTMcenterState state)
{
  const int niter_tot = niter*thin + nburn;
  int nsaved{0};
  sample::GSL_RNG engine(param.seed);

  check_positive_finite_centers(param.gamma, "gamma");
  check_positive_finite_centers(param.beta, "beta");
  if(param.sigma <= 0.0 || param.sigma >= 1.0 || !std::isfinite(param.sigma))
    throw std::runtime_error("Error in GibbsSampler_DTM_centers_c_core: sigma must be in (0,1)");
  check_center_state_shapes(state, H, V, Ttot);

  Rcpp::List Xi_mcmc(niter);
  Rcpp::List S_mcmc(niter);
  Rcpp::List U_mcmc(niter);
  Rcpp::List N_mcmc(niter);
  Rcpp::List Lambda_star_mcmc(niter);
  Rcpp::List Lambda_star_by_center_mcmc(niter);
  Rcpp::List Zeta_mcmc(niter);
  Rcpp::List centers_aux_mcmc(niter);
  VecUnsCol M_mcmc{VecUnsCol::Zero(niter)};
  VecUnsCol Mstar_mcmc{VecUnsCol::Zero(niter)};
  VecCol omega_mcmc{VecCol::Zero(niter)};
  VecCol t_sigma_gamma_mcmc{VecCol::Zero(niter)};

  Progress progress_bar(niter_tot, param.print);
  for(int it = 0; it < niter_tot; it++){
    const double t_sigma_gamma = compute_t_sigma_gamma_centers(param.gamma, param.sigma, H, state.M);

    for(unsigned int m = 0; m < state.M; m++){
      if(param.UpdateXi)
        state.Xi[m] = sample_Xi_tl_centers(engine, state.Xi[m], state.S[m], state.N[m],
                                           PROCESS_PHI, param.sigma, param.beta, t_sigma_gamma);

      if(param.UpdateS)
        state.S[m] = sample_Stl(engine, state.Xi[m], state.U[m],
                                PROCESS_PHI, param.sigma, param.beta);

      if(param.UpdateU)
        state.U[m] = sample_Utl(engine, state.S[m], t_sigma_gamma);

      if(param.UpdateLambda)
        state.Lambda[m] = sample_Lambda_itlm(engine, state.Dl[m], state.Xi[m], state.Zeta[m]);
    }

    if(param.UpdateDitl){
      auto aux_D = sample_Ditlm(engine, state.Lambda, state.Xi, D);
      state.Dl = aux_D.first;
      state.N = aux_D.second;
    }

    DTMcenterLambdaStar Lambda_star = build_Lambda_star_centers(state.Lambda, state.Xi);
    ClusTopicUpdate centers_update;
    if(param.UpdateCenters){
      ClusTopicParams clus_param = param.clus;
      clus_param.omega = state.omega;
      centers_update = sample_ClusTopic_partition_fixedM(engine, Lambda_star.all,
                                                         allocated_zeta(state.Zeta, state.M),
                                                         clus_param);
      state.Zeta = centers_update.Zeta;
      state.Mstar = 0;
      state.omega = centers_update.omega;

      if(param.UpdateDitl){
        auto aux_D = sample_Ditlm(engine, state.Lambda, state.Xi, D);
        state.Dl = aux_D.first;
        state.N = aux_D.second;
      }
    }

    if(it >= nburn && (it - nburn)%thin == 0){
      if(nsaved >= niter)
        throw std::runtime_error("Error in GibbsSampler_DTM_centers_c_core: too many saved states");

      Xi_mcmc[nsaved] = Rcpp::wrap(state.Xi);
      S_mcmc[nsaved] = Rcpp::wrap(state.S);
      U_mcmc[nsaved] = Rcpp::wrap(state.U);
      N_mcmc[nsaved] = Rcpp::wrap(state.N);
      Lambda_star_mcmc[nsaved] = Lambda_star.all;
      Lambda_star_by_center_mcmc[nsaved] = Rcpp::wrap(Lambda_star.by_center);
      Zeta_mcmc[nsaved] = Rcpp::wrap(state.Zeta);

      if(param.UpdateCenters){
        centers_aux_mcmc[nsaved] = Rcpp::List::create(
          Rcpp::Named("c") = centers_update.c,
          Rcpp::Named("c_raw") = centers_update.aux.c_raw,
          Rcpp::Named("allocation_prob") = centers_update.aux.allocation_prob,
          Rcpp::Named("cluster_size") = centers_update.aux.cluster_size,
          Rcpp::Named("A") = centers_update.aux.A,
          Rcpp::Named("mstar_prob") = centers_update.aux.mstar_prob,
          Rcpp::Named("log_mstar_prob") = centers_update.aux.log_mstar_prob,
          Rcpp::Named("phi") = centers_update.aux.phi,
          Rcpp::Named("delta") = Rcpp::wrap(centers_update.aux.delta),
          Rcpp::Named("log_acc_zeta") = centers_update.aux.log_acc_zeta,
          Rcpp::Named("accept_zeta") = centers_update.aux.accept_zeta,
          Rcpp::Named("omega_shape") = centers_update.aux.omega_shape,
          Rcpp::Named("omega_rate") = centers_update.aux.omega_rate
        );
      }
      else{
        centers_aux_mcmc[nsaved] = Rcpp::List::create();
      }

      M_mcmc(nsaved) = state.M;
      Mstar_mcmc(nsaved) = state.Mstar;
      omega_mcmc(nsaved) = state.omega;
      t_sigma_gamma_mcmc(nsaved) = compute_t_sigma_gamma_centers(param.gamma, param.sigma, H, state.M);
      nsaved++;
    }

    try{
      Rcpp::checkUserInterrupt();
    }
    catch(Rcpp::internal::InterruptedException e){
      throw std::runtime_error("Execution stopped by the user");
    }
    progress_bar.increment();
  }

  return Rcpp::List::create(
    Rcpp::Named("Xi") = Xi_mcmc,
    Rcpp::Named("S") = S_mcmc,
    Rcpp::Named("U") = U_mcmc,
    Rcpp::Named("N") = N_mcmc,
    Rcpp::Named("Lambda_star") = Lambda_star_mcmc,
    Rcpp::Named("Lambda_star_by_center") = Lambda_star_by_center_mcmc,
    Rcpp::Named("Zeta") = Zeta_mcmc,
    Rcpp::Named("M") = M_mcmc,
    Rcpp::Named("Mstar") = Mstar_mcmc,
    Rcpp::Named("omega") = omega_mcmc,
    Rcpp::Named("gamma") = param.gamma,
    Rcpp::Named("sigma") = param.sigma,
    Rcpp::Named("beta") = param.beta,
    Rcpp::Named("phi_process") = PROCESS_PHI,
    Rcpp::Named("phi_centers") = param.clus.phi,
    Rcpp::Named("t_sigma_gamma") = t_sigma_gamma_mcmc,
    Rcpp::Named("centers_aux") = centers_aux_mcmc
  );
}
