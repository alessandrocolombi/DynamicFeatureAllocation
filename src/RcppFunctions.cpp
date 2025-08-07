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
//	log_stable_sum
//------------------------------------------------------------------------------------------------------------------------------------------------------

double log_stable_sum(const Rcpp::NumericVector& a, const bool is_log, const double& val_max){

	double inf = std::numeric_limits<double>::infinity();

	if(a.size() == 0)
		return 0.0;

	// Do not checks if it is really the max value
	if(is_log){ // a contains values in log scale

		if(val_max == -inf) // need to handle the case where all values are -inf
			return -inf;

		return (val_max +
					std::log( std::accumulate(   a.cbegin(), a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )   )
			   );
	}
	else{

		if(val_max < 0)
			throw std::runtime_error("log_stable_sum, if is_log is FALSE, the maximum value can not be negative ");
		if(val_max == 0)
			return 0.0;

		return ( std::log(val_max) +
				 std::log( std::accumulate(   a.cbegin(), a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   ) )
			   );
	}

}

// In this version of the formula, the maximum value is computed
// [[Rcpp::export]]
double log_stable_sum(const Rcpp::NumericVector& a, const bool is_log){
	if(a.size() == 0)
		return 0.0;

	// Computes maximum value
	auto it_max{std::max_element(a.cbegin(), a.cend())};
	double val_max{*it_max};
	// Calls the specialized version
	return log_stable_sum(a,is_log,val_max);
}


// [[Rcpp::export]]
MyTraits::MatCol sample_A( const int& K, const MyTraits::VecCol& x, const MyTraits::VecCol& mu0, 
													 const double& sig2_X, const double& sig2_A, const int& seed  )
{

	int D(x.size() ); // problem size
	if(K <= 0)
		return(MyTraits::MatCol(0,D));
	if(D <= 0)
		throw std::runtime_error("Error in sample_A: invalid number of cols (D)");

	if(seed == 0)
		throw std::runtime_error("Error in sample_A: the seed can not be equal to 0");


	// Define basic quantities
	sample::GSL_RNG random_engine(seed); // GSL random engine to sample from random distribution
	MyTraits::MatCol Ones_mat(MyTraits::MatCol::Constant(K,K,1.0));
	MyTraits::VecCol Ones_vec(MyTraits::VecCol::Constant(K,1.0));
	MyTraits::MatCol Id(MyTraits::MatCol::Identity(K,K));
	double sig_ratio = sig2_X/sig2_A;

	// Compute posterior quantities
	MyTraits::MatCol Omega(Ones_mat + sig_ratio*Id); // Precision matrix
	MyTraits::MatCol Omega_inv( 1.0/sig_ratio*(Id - sig2_A/(sig2_X + (double)K*sig2_A)*Ones_mat) ); // Covariance matrix
	//MyTraits::MatCol Atilde = Omega_inv * Ones_vec * (x + sig_ratio*mu0).transpose() ; // mean

	// The mean of A|X is Atilde, that is a (KxD) matrix such that
	// A[j,] = 1/(sig2_X + K*sig2_A) * (sig2_A*x + sig2_X*mu0), for j = 1,...,K 
	MyTraits::VecRow Atilde_j = 1/(sig2_X + (double)K*sig2_A) * (sig2_A*x + sig2_X*mu0);

	// Sample
	sample::rmvnorm rmv; //Covariance parametrization
	MyTraits::MatCol Anew(MyTraits::MatCol::Constant(K,D,0.0));
	for(std::size_t ii=0; ii < D; ii++ ){
		MyTraits::VecCol mean = Atilde_j(ii) * Ones_vec; // Atilde_j(ii) is repeated K times
		MyTraits::MatCol Sigma{sig2_X*Omega_inv}; // scale covariance matrix
		Anew.col(ii) = rmv(random_engine, mean, Sigma); // j-th element for all K features
	}

	return( Anew );
}

// [[Rcpp::export]]
double log_dmarg_img( const int& K, const MyTraits::VecCol& x, const MyTraits::VecCol& mu0, 
											const double& sig2_X, const double& sig2_A)
{
	if(sig2_X <= 0)
		throw std::runtime_error("Error in log_dmarg_img: Negative sig2_X");
	if(sig2_A <= 0)
		throw std::runtime_error("Error in log_dmarg_img: Negative sig2_A");

	int D(x.size()); // problem size

	MyTraits::VecCol mean{MyTraits::VecCol::Constant(D,0.0)};
	double scale = sig2_X;
	if(K > 0){
		mean += (double)K * mu0;
		scale += sig2_A*(double)K;
	}
	if(scale <= 0)
		throw std::runtime_error("Error in log_dmarg_img: Negative variance");

	double res = 0;
	for(std::size_t j=0; j < D; j++ ){
		res += log_dnorm(x[j], mean[j], std::sqrt(scale) );
	}
	return res;
}

// --------------------------------------------------------------------------------------------
// Test functions
// --------------------------------------------------------------------------------------------


Rcpp::NumericVector prova(Rcpp::NumericVector x)
{
  return x+x;
}

