#ifndef __FULLCONDITIONALS_HPP__
#define __FULLCONDITIONALS_HPP__

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>

// Include file with basic libraries to include
#include "headers.h"
#include "recurrent_traits.h"

#include "utils.h"
#include "mysample.h"


using namespace Rcpp;

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Sampling
//------------------------------------------------------------------------------------------------------------------------------------------------------

int sample_Ditl(sample::GSL_RNG const & engine, const std::vector<MatCol>& Lambda_itl, const MatIntCol& Xi, const MatIntCol& D);

int sample_Lambda_itl(sample::GSL_RNG const & engine, const std::vector<MatCol>& D_itl, const MatIntCol& Xi, const double& delta);

int sample_Stl(sample::GSL_RNG const & engine, const MatIntCol& Xi, const MatCol& U, const double& phi, const double& sigma, const double& b);

int sample_Utl(sample::GSL_RNG const & engine, const MatCol& S, const double& t_sigma_gamma);

int sample_Xi_tl(sample::GSL_RNG const & engine, const MatIntCol& Xi_old, const MatCol& S, const MatUnsCol& N_tl, 
                 const double& phi, const double& sigma, const double& b, const double& t_sigma_gamma);

int sample_hyparams(sample::GSL_RNG const & engine, const MatIntCol& Xi, const MatCol& S, 
                    const double& phi_old, const double& gamma_old, const double& sigma_old, const double& b_old, 
                    const double& t_sigma_gamma_old,
                    const double& a_phi, const double& b_phi, 
                    const double& a_gamma, const double& b_gamma, 
                    const double& a_sigma, const double& b_sigma, 
                    const double& a_beta, const double& b_beta, 
                    const double& var_phi, const double& var_gamma, const double& var_sigma, const double& var_beta,
                    bool UpdatePhi, bool UpdateGamma, bool UpdateSigma, bool UpdateBeta);
//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Utilities
//------------------------------------------------------------------------------------------------------------------------------------------------------
void find_indices(const VecIntCol& Z, std::vector<int>& idx_born, std::vector<int>& idx_surv, std::vector<int>& idx_noact);

#endif
