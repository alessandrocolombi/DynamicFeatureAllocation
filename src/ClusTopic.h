#ifndef __CLUSTOPIC_HPP__
#define __CLUSTOPIC_HPP__

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>

#include "headers.h"
#include "recurrent_traits.h"
#include "mysample.h"

using namespace Rcpp;

//------------------------------------------------------------------------------------------------------------------------------------------------------
//  Data structures for random-center topic clustering
//------------------------------------------------------------------------------------------------------------------------------------------------------

struct ClusTopicParams
{
  double phi{1.0};     // fixed common concentration: zeta_m = phi * delta_m
  double delta0{1.0};  // symmetric Dirichlet hyperparam. for delta_m
  double omega{1.0};   // hyparam. for Poisson distributed number of clusters
  double a_omega{1.0}; // shape hyparam. for omega
  double b_omega{1.0}; // rate hyparam. for omega
  double var_delta{0.01}; // adaptive variance in proposal for delta
  unsigned int mstar_max{100}; // truncation level
  bool UpdateZeta{true};
  bool UpdateOmega{false};
};

struct ClusTopicZetaDraw
{
  VecCol zeta;
  double phi{1.0};
  VecCol delta;
};

struct ClusTopicZetaMH
{
  VecCol zeta;
  double phi{1.0};
  VecCol delta;
  double log_acc{0.0};
  bool accepted{false};
};

struct ClusTopicAux
{
  VecIntCol c_raw;
  MatCol allocation_prob;
  VecUnsCol cluster_size;
  MatCol A;
  VecCol mstar_prob;
  VecCol log_mstar_prob;
  double phi{1.0};
  std::vector<VecCol> delta;
  VecCol log_acc_zeta;
  VecUnsCol accept_zeta;
  double omega_shape{1.0};
  double omega_rate{1.0};
};

struct ClusTopicUpdate
{
  VecIntCol c;
  unsigned int M{0};
  unsigned int Mstar{0};
  double omega{1.0};
  std::vector<VecCol> Zeta;
  ClusTopicAux aux;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------
//  Utilities and sampling steps
//------------------------------------------------------------------------------------------------------------------------------------------------------

double log_ClusTopic_dirichlet_density(const VecCol& lambda, const VecCol& zeta);

double log_ClusTopic_allocated_zeta_full_conditional(const VecCol& A_m, const unsigned int& n_m,
                                                     const double& phi, const VecCol& delta,
                                                     const double& delta0);

VecUnsCol compute_ClusTopic_cluster_sizes(const VecIntCol& c, const unsigned int& M);

MatCol compute_ClusTopic_A(const MatCol& Lambda_star, const VecIntCol& c, const unsigned int& M);

ClusTopicZetaDraw split_ClusTopic_zeta(const VecCol& zeta, const double& phi);

ClusTopicZetaDraw sample_ClusTopic_prior_zeta(sample::GSL_RNG const & engine, const unsigned int& V,
                                              const double& phi, const double& delta0);

unsigned int sample_ClusTopic_Mstar(sample::GSL_RNG const & engine, const unsigned int& M,
                                    const unsigned int& K, const double& omega,
                                    const unsigned int& mstar_max,
                                    VecCol& mstar_prob, VecCol& log_mstar_prob);

ClusTopicZetaMH sample_ClusTopic_allocated_zeta(sample::GSL_RNG const & engine, const VecCol& zeta_old,
                                                const VecCol& A_m, const unsigned int& n_m,
                                                const double& phi, const double& delta0,
                                                const double& var_delta);

ClusTopicUpdate sample_ClusTopic_partition(sample::GSL_RNG const & engine, const MatCol& Lambda_star,
                                           const std::vector<VecCol>& Zeta_old,
                                           const ClusTopicParams& param);

ClusTopicUpdate sample_ClusTopic_partition_fixedM(sample::GSL_RNG const & engine, const MatCol& Lambda_star,
                                                  const std::vector<VecCol>& Zeta_old,
                                                  const ClusTopicParams& param);

#endif
