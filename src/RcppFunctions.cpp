// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]

// Include file with basic libraries to include
#include "headers.h"
#include "recurrent_traits.h"

#include "utils.h"
#include "mysample.h"

#include "FC.h"


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
													 const double& sig2_X, const double& sig2_A, const unsigned int& seed  )
{
	
	int D(x.size() ); // problem size
	if(K <= 0)
		return(MyTraits::MatCol(0,D));
	if(D <= 0)
		throw std::runtime_error("Error in sample_A: invalid number of cols (D)");

	if(seed <= 0)
		throw std::runtime_error("Error in sample_A: the seed must be strictly positive");


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
// Truncated Gibbs Sampling for Topic Modeling
// --------------------------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List GibbsSampler_DTM_c(const int& niter, const int& nburn, 
															const MatIntCol& D, const int& H, const int& V, const int& Ttot,
															const Rcpp::List& param_DTM, const Rcpp::List& init_DTM)
{
	const int niter_tot = niter + nburn;
	//Read param_DTM
	double delta     = as<double>(param_DTM["delta"]);
	double a_phi     = as<double>(param_DTM["a_phi"]);
	double b_phi     = as<double>(param_DTM["b_phi"]);
	double a_gamma   = as<double>(param_DTM["a_gamma"]);
	double b_gamma   = as<double>(param_DTM["b_gamma"]);
	double a_sigma   = as<double>(param_DTM["a_sigma"]);
	double b_sigma   = as<double>(param_DTM["b_sigma"]);
	double a_beta    = as<double>(param_DTM["a_beta"]);
	double b_beta    = as<double>(param_DTM["b_beta"]);
	
	double var_phi   = as<double>(param_DTM["var_phi"]);
	double var_gamma = as<double>(param_DTM["var_gamma"]);
	double var_sigma = as<double>(param_DTM["var_sigma"]);
	double var_beta  = as<double>(param_DTM["var_beta"]);
	
	bool UpdateDitl   = as<bool>(param_DTM["UpdateDitl"]);
	bool UpdateS      = as<bool>(param_DTM["UpdateS"]);
	bool UpdateLambda = as<bool>(param_DTM["UpdateLambda"]);
	bool UpdateXi     = as<bool>(param_DTM["UpdateXi"]);
	bool UpdateU      = as<bool>(param_DTM["UpdateU"]);
	bool UpdatePhi    = as<bool>(param_DTM["UpdatePhi"]);
	bool UpdateGamma  = as<bool>(param_DTM["UpdateGamma"]);
	bool UpdateSigma  = as<bool>(param_DTM["UpdateSigma"]);
	bool UpdateBeta   = as<bool>(param_DTM["UpdateBeta"]);
	bool print        = as<bool>(param_DTM["print"]);
	
	int seed   = as<int>(param_DTM["seed"]);
	sample::GSL_RNG engine(seed);

	// Read initial values
  MatIntCol Xi0 = init_DTM["Xi0"];
  MatCol S0     = init_DTM["S0"];
  List Lambda0_list = init_DTM["Lambda0"];
  std::vector<MatCol> Lambda0(H);
  for (int l = 0; l < H; l++) {
    Lambda0[l] = as<MatCol>(Lambda0_list[l]);
  }

  // Check Lambda0
  for(int l=0; l < H; l++){
  	VecCol aux = Lambda0[l].colwise().sum();
  	if (!((aux.array() - 1.0).abs() <= 1e-10).all()) {
  	    throw std::runtime_error("Error in Lambda0: column values must sum to 1");
  	}
  }

  double phi0   = as<double>(init_DTM["phi0"]);
  double gamma0 = as<double>(init_DTM["gamma0"]);
  double sigma0 = as<double>(init_DTM["sigma0"]);
  double beta0  = as<double>(init_DTM["beta0"]);
  double t_sigma_gamma0 = std::exp( 1.0/sigma0 * ( std::log(sigma0*(double)H) - std::log(gamma0) )  );

  // Initilize Dl and U
  MatCol U0 = sample_Utl(engine, S0, t_sigma_gamma0);
  auto temp = sample_Ditl(engine, Lambda0, Xi0, D);
  std::vector<MatUnsCol> Dl0 = temp.first;
  MatUnsCol N0 = temp.second;

  // Main objects initialization
  std::vector<MatIntCol> Xi_mcmc(niter_tot+1, MatIntCol::Zero(H,Ttot)); Xi_mcmc[0] = Xi0;
  std::vector<MatCol> S_mcmc(niter_tot+1, MatCol::Zero(H,Ttot));        S_mcmc[0] = S0;
  std::vector<MatCol> U_mcmc(niter_tot+1,  MatCol::Zero(H,Ttot));       U_mcmc[0] = U0;
  std::vector<MatUnsCol> N_mcmc(niter_tot+1,  MatUnsCol::Zero(Ttot,H)); N_mcmc[0] = N0;
  std::vector<std::vector<MatUnsCol>> Dl_mcmc(niter_tot+1, Dl0);
  std::vector<std::vector<MatCol>> Lambda_mcmc(niter_tot+1, Lambda0);

  std::vector<double> phi_mcmc(niter_tot+1,-1.0);   phi_mcmc[0]   = phi0;
  std::vector<double> gamma_mcmc(niter_tot+1,-1.0); gamma_mcmc[0] = gamma0;
  std::vector<double> sigma_mcmc(niter_tot+1,-1.0); sigma_mcmc[0] = sigma0;
  std::vector<double> beta_mcmc(niter_tot+1,-1.0);  beta_mcmc[0] = beta0;
  std::vector<double> t_sigma_gamma_mcmc(niter_tot+1,-1.0);  t_sigma_gamma_mcmc[0] = t_sigma_gamma0;

  // Start MCMC loop
  Rcpp::Rcout<<"Preprocessing finished. Start MCMC ... "<<std::endl;
  Progress progress_bar(niter_tot, print); // Initialize progress bar
  for(int it = 1; it <= niter_tot; it++){

  	// ----------------------------------
  	if(UpdateXi){
  		Xi_mcmc[it] = sample_Xi_tl(engine, Xi_mcmc[it-1], S_mcmc[it-1], N_mcmc[it-1],
  		                       			phi_mcmc[it-1], sigma_mcmc[it-1], beta_mcmc[it-1], 
  		                       			t_sigma_gamma_mcmc[it-1]);
  	}
  	else{
  		Xi_mcmc[it] = Xi_mcmc[it-1];
  	}
		// ----------------------------------
  	if(UpdateS){
			S_mcmc[it] = sample_Stl(engine, Xi_mcmc[it], U_mcmc[it-1], phi_mcmc[it-1], sigma_mcmc[it-1], beta_mcmc[it-1]);
  	}
  	else{
  		S_mcmc[it] = S_mcmc[it-1];
  	}
  	// ----------------------------------
  	if(UpdateU){
			U_mcmc[it] = sample_Utl(engine, S_mcmc[it], t_sigma_gamma_mcmc[it-1]);
  	}
  	else{
  		U_mcmc[it] = U_mcmc[it-1];
  	}
  	// ----------------------------------
  	if(UpdateLambda){
  		Lambda_mcmc[it-1] = sample_Lambda_itl(engine, Dl_mcmc[it-1], Xi_mcmc[it], delta);
  	}
  	else{
  		Lambda_mcmc[it] = Lambda_mcmc[it-1];
  	}
  	// Check Lambda_mcmc[it]
  	// Is this necessary??
  	for(int l=0; l < H; l++){
  		VecCol aux = Lambda_mcmc[it][l].colwise().sum();
  		if (!((aux.array() - 1.0).abs() <= 1e-10).all()) {
  		    throw std::runtime_error("Error in Lambda: column values must sum to 1");
  		}
  	}
  	// ----------------------------------
  	if(UpdateDitl){
  		auto aux = sample_Ditl(engine, Lambda_mcmc[it], Xi_mcmc[it], D);
			Dl_mcmc[it] = aux.first;
			N_mcmc[it]  = aux.second;
  	}
  	else{
  		Dl_mcmc[it] = Dl_mcmc[it-1];
  		N_mcmc[it]  = N_mcmc[it-1];
  	}
  	// ----------------------------------
  	if(UpdatePhi || UpdateBeta || UpdateGamma || UpdateSigma){
  		VecCol aux = sample_hyparams( engine, Xi_mcmc[it], S_mcmc[it], 
                         						phi_mcmc[it-1],  gamma_mcmc[it-1],  sigma_mcmc[it-1],  beta_mcmc[it-1], t_sigma_gamma_mcmc[it-1],
                         						a_phi, b_phi, a_gamma, b_gamma, a_sigma, b_sigma, a_beta, b_beta,  
                         						var_phi,  var_gamma,  var_sigma,  var_beta,
			                        			UpdatePhi, UpdateGamma, UpdateSigma, UpdateBeta);
	  	phi_mcmc[it]   = aux[0];
	  	gamma_mcmc[it] = aux[1];
	  	sigma_mcmc[it] = aux[2];
	  	beta_mcmc[it]  = aux[3]; 
	  	t_sigma_gamma_mcmc[it] = aux[4];
  	}
  	else{
  		phi_mcmc[it]   = phi_mcmc[it-1];
  		gamma_mcmc[it] = gamma_mcmc[it-1];
  		sigma_mcmc[it] = sigma_mcmc[it-1];
  		beta_mcmc[it]  = beta_mcmc[it-1]; 
  		t_sigma_gamma_mcmc[it] = t_sigma_gamma_mcmc[it-1];
  	}

  	//Check for User Interruption
    try{
    	Rcpp::checkUserInterrupt();
    }
    catch(Rcpp::internal::InterruptedException e){ 
    	//Print error and return
      throw std::runtime_error("Execution stopped by the user");
    }
  	progress_bar.increment(); //update progress bar
  }

  return Rcpp::List::create(  
  	Rcpp::Named("Xi") = Xi_mcmc,
  	Rcpp::Named("S") = S_mcmc,
  	Rcpp::Named("U") = U_mcmc,
  	Rcpp::Named("N") = N_mcmc,
  	Rcpp::Named("Dl") = Dl_mcmc,
  	Rcpp::Named("Lambda") = Lambda_mcmc,
  	Rcpp::Named("phi") = phi_mcmc,
  	Rcpp::Named("gamma") = gamma_mcmc,
  	Rcpp::Named("sigma") = sigma_mcmc,
  	Rcpp::Named("beta") = beta_mcmc,
  	Rcpp::Named("t_sigma_gamma") = t_sigma_gamma_mcmc
  );
}

// --------------------------------------------------------------------------------------------
// Test functions
// --------------------------------------------------------------------------------------------


// [[Rcpp::export]]
MatIntCol sample_Xi_tl(const int& seed, const MatIntCol& Xi_old, const MatCol& S, const MatUnsCol& N_tl, 
                       const double& phi, const double& sigma, const double& b, const double& t_sigma_gamma)
{
	sample::GSL_RNG engine(seed);
	return sample_Xi_tl(engine, Xi_old, S, N_tl, phi, sigma, b, t_sigma_gamma);
}

Rcpp::NumericVector prova(Rcpp::NumericVector x)
{
  return x+x;
}


