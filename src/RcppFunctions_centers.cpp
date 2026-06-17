// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>
#include <progress.hpp>
#include <progress_bar.hpp>

#include "headers.h"
#include "recurrent_traits.h"
#include "utils.h"
#include "mysample.h"
#include "FC.h"
#include "ClusTopic.h"
#include "DTMcenters.h"

using namespace Rcpp;

namespace {

double get_double_or_default_centers(const Rcpp::List& x, const char* name, const double& default_value)
{
  if(x.containsElementNamed(name))
    return Rcpp::as<double>(x[name]);
  return default_value;
}

int get_int_or_default_centers(const Rcpp::List& x, const char* name, const int& default_value)
{
  if(x.containsElementNamed(name))
    return Rcpp::as<int>(x[name]);
  return default_value;
}

bool get_bool_or_default_centers(const Rcpp::List& x, const char* name, const bool& default_value)
{
  if(x.containsElementNamed(name))
    return Rcpp::as<bool>(x[name]);
  return default_value;
}

ClusTopicParams read_ClusTopic_params_centers(const Rcpp::List& param)
{
  ClusTopicParams clus;
  clus.a_phi       = get_double_or_default_centers(param, "a_phi_centers", clus.a_phi);
  clus.b_phi       = get_double_or_default_centers(param, "b_phi_centers", clus.b_phi);
  clus.delta0      = get_double_or_default_centers(param, "delta0_centers", clus.delta0);
  clus.omega       = get_double_or_default_centers(param, "omega", clus.omega);
  clus.a_omega     = get_double_or_default_centers(param, "a_omega", clus.a_omega);
  clus.b_omega     = get_double_or_default_centers(param, "b_omega", clus.b_omega);
  clus.var_phi     = get_double_or_default_centers(param, "var_phi_centers", clus.var_phi);
  clus.var_delta   = get_double_or_default_centers(param, "var_delta_centers", clus.var_delta);
  clus.mstar_max   = (unsigned int)get_int_or_default_centers(param, "mstar_max", clus.mstar_max);
  clus.UpdateZeta  = get_bool_or_default_centers(param, "UpdateCenters", clus.UpdateZeta);
  clus.UpdateOmega = get_bool_or_default_centers(param, "UpdateOmega", clus.UpdateOmega);
  return clus;
}

DTMcenterParams read_DTMcenter_params(const Rcpp::List& param)
{
  DTMcenterParams out;
  out.gamma        = get_double_or_default_centers(param, "gamma", out.gamma);
  out.sigma        = get_double_or_default_centers(param, "sigma", out.sigma);
  out.beta         = get_double_or_default_centers(param, "beta", out.beta);
  out.UpdateDitl   = get_bool_or_default_centers(param, "UpdateDitl", out.UpdateDitl);
  out.UpdateS      = get_bool_or_default_centers(param, "UpdateS", out.UpdateS);
  out.UpdateLambda = get_bool_or_default_centers(param, "UpdateLambda", out.UpdateLambda);
  out.UpdateXi     = get_bool_or_default_centers(param, "UpdateXi", out.UpdateXi);
  out.UpdateU      = get_bool_or_default_centers(param, "UpdateU", out.UpdateU);
  out.UpdateCenters = get_bool_or_default_centers(param, "UpdateCenters", out.UpdateCenters);
  out.print        = get_bool_or_default_centers(param, "print", out.print);
  out.seed         = (unsigned int)get_int_or_default_centers(param, "seed", out.seed);
  out.clus         = read_ClusTopic_params_centers(param);
  return out;
}

std::vector<MatIntCol> read_MatIntCol_list(const Rcpp::List& x)
{
  std::vector<MatIntCol> out(x.size());
  for(int m = 0; m < x.size(); m++)
    out[m] = Rcpp::as<MatIntCol>(x[m]);
  return out;
}

std::vector<MatCol> read_MatCol_list(const Rcpp::List& x)
{
  std::vector<MatCol> out(x.size());
  for(int m = 0; m < x.size(); m++)
    out[m] = Rcpp::as<MatCol>(x[m]);
  return out;
}

std::vector<VecCol> read_VecCol_list(const Rcpp::List& x)
{
  std::vector<VecCol> out(x.size());
  for(int m = 0; m < x.size(); m++)
    out[m] = Rcpp::as<VecCol>(x[m]);
  return out;
}

std::vector<std::vector<MatCol>> read_nested_Lambda_list(const Rcpp::List& x)
{
  std::vector<std::vector<MatCol>> out(x.size());
  for(int m = 0; m < x.size(); m++){
    Rcpp::List x_m = x[m];
    out[m].resize(x_m.size());
    for(int l = 0; l < x_m.size(); l++)
      out[m][l] = Rcpp::as<MatCol>(x_m[l]);
  }
  return out;
}

DTMcenterState read_DTMcenter_init(const Rcpp::List& init, const MatIntCol& D,
                                   const DTMcenterParams& param,
                                   const int& H, const int& V, const int& Ttot)
{
  DTMcenterState state;
  state.Xi = read_MatIntCol_list(init["Xi0"]);
  state.S = read_MatCol_list(init["S0"]);
  state.Lambda = read_nested_Lambda_list(init["Lambda0"]);
  state.Zeta = read_VecCol_list(init["Zeta0"]);
  state.M = state.Xi.size();
  state.Mstar = 0;
  state.omega = param.clus.omega;

  if(state.Zeta.size() != state.M)
    throw std::runtime_error("Error in read_DTMcenter_init: fixed-M sampler requires length(Zeta0) == M0");

  if(init.containsElementNamed("U0")){
    state.U = read_MatCol_list(init["U0"]);
  }
  else{
    state.U.resize(state.M);
    const double t_sigma_gamma = compute_t_sigma_gamma_centers(param.gamma, param.sigma, H, state.M);
    for(unsigned int m = 0; m < state.M; m++)
      state.U[m] = sample_Utl(sample::GSL_RNG(param.seed + m + 1), state.S[m], t_sigma_gamma);
  }

  if(init.containsElementNamed("omega0"))
    state.omega = Rcpp::as<double>(init["omega0"]);

  auto aux_D = sample_Ditlm(sample::GSL_RNG(1), state.Lambda, state.Xi, D);
  state.Dl = aux_D.first;
  state.N = aux_D.second;
  return state;
}

} // anonymous namespace


// [[Rcpp::export]]
Rcpp::List GibbsSampler_DTM_centers_c(const int& niter, const int& nburn, const int& thin,
                                      const MatIntCol& D, const int& H, const int& V,
                                      const int& Ttot, const Rcpp::List& param_DTM_centers,
                                      const Rcpp::List& init_DTM_centers)
{
  DTMcenterParams param = read_DTMcenter_params(param_DTM_centers);
  DTMcenterState state = read_DTMcenter_init(init_DTM_centers, D, param, H, V, Ttot);
  return GibbsSampler_DTM_centers_c_core(niter, nburn, thin, D, H, V, Ttot, param, state);
}
