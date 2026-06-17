#ifndef __DTMCENTERS_HPP__
#define __DTMCENTERS_HPP__

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
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

using namespace Rcpp;

struct DTMcenterParams
{
  double gamma{1.0};
  double sigma{0.5};
  double beta{1.0};
  bool UpdateDitl{true};
  bool UpdateS{true};
  bool UpdateLambda{true};
  bool UpdateXi{true};
  bool UpdateU{true};
  bool UpdateCenters{true};
  bool print{true};
  unsigned int seed{1};
  ClusTopicParams clus;
};

struct DTMcenterState
{
  std::vector<MatIntCol> Xi;
  std::vector<MatCol> S;
  std::vector<MatCol> U;
  std::vector<std::vector<MatCol>> Lambda;
  std::vector<std::vector<MatUnsCol>> Dl;
  std::vector<MatUnsCol> N;
  std::vector<VecCol> Zeta;
  unsigned int M{0};
  unsigned int Mstar{0};
  double omega{1.0};
};

struct DTMcenterLambdaStar
{
  std::vector<MatCol> by_center;
  MatCol all;
  VecIntCol center_id;
  VecIntCol row_id;
  VecIntCol born_time;
};

double compute_t_sigma_gamma_centers(const double& gamma, const double& sigma,
                                     const unsigned int& H, const unsigned int& M);

MatIntCol sample_Xi_tl_centers(sample::GSL_RNG const & engine,
                               const MatIntCol& Xi_old,
                               const MatCol& S,
                               const MatUnsCol& N_tl,
                               const double& phi,
                               const double& sigma,
                               const double& b,
                               const double& t_sigma_gamma);

std::vector<MatCol> sample_Lambda_itlm(sample::GSL_RNG const & engine,
                                       const std::vector<MatUnsCol>& D_itl,
                                       const MatIntCol& Xi,
                                       const VecCol& zeta);

std::pair<std::vector<std::vector<MatUnsCol>>, std::vector<MatUnsCol>>
sample_Ditlm(sample::GSL_RNG const & engine,
             const std::vector<std::vector<MatCol>>& Lambda,
             const std::vector<MatIntCol>& Xi,
             const MatIntCol& D);

DTMcenterLambdaStar build_Lambda_star_centers(const std::vector<std::vector<MatCol>>& Lambda,
                                              const std::vector<MatIntCol>& Xi);

void resize_DTMcenter_state_after_clustering(sample::GSL_RNG const & engine,
                                             DTMcenterState& state,
                                             const unsigned int& M_new,
                                             const unsigned int& H,
                                             const unsigned int& V,
                                             const unsigned int& Ttot,
                                             const double& gamma,
                                             const double& sigma,
                                             const double& beta);

Rcpp::List GibbsSampler_DTM_centers_c_core(const int& niter, const int& nburn, const int& thin,
                                           const MatIntCol& D, const int& H, const int& V,
                                           const int& Ttot, const DTMcenterParams& param,
                                           DTMcenterState state);

#endif
