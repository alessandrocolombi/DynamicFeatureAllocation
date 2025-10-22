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

int sample_Lambda_itl(sample::GSL_RNG const & engine, const std::vector<MatCol>& D_itl, const MatIntCol& Xi);

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Utilities
//------------------------------------------------------------------------------------------------------------------------------------------------------
void find_indices(const VecIntCol& Z, std::vector<int>& idx_born, std::vector<int>& idx_surv, std::vector<int>& idx_noact);

#endif
