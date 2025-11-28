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
Rcpp::List GibbsSampler_DTM_c(const int& niter, const int& nburn, const int& thin,
															const MatIntCol& D, const int& H, const int& V, const int& Ttot,
															const Rcpp::List& param_DTM, const Rcpp::List& init_DTM)
{
	const int niter_tot = niter*thin + nburn;
	int nsaved{0};

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
	bool JointAdp     = as<bool>(param_DTM["JointAdp"]);
	
	int seed   = as<int>(param_DTM["seed"]);
	sample::GSL_RNG engine(seed);

	// Read initial values
	Rcpp::Rcout<<"Start Preprocessing: read initial values ... ";
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
  Rcpp::Rcout<<" compute U0 and D_itl0 ... ";
  MatCol U0 = sample_Utl(engine, S0, t_sigma_gamma0);
  auto temp = sample_Ditl(engine, Lambda0, Xi0, D);
  std::vector<MatUnsCol> Dl0 = temp.first;
  MatUnsCol N0 = temp.second;

  Rcpp::Rcout<<" main objects definition ... ";
   // Main objects initialization
  MatIntCol Xi = Xi0;
  MatCol S = S0;
  MatCol U = U0;
  MatUnsCol N = N0;
  std::vector<MatUnsCol> Dl = Dl0;
  std::vector<MatCol> Lambda = Lambda0;

  double phi = phi0;
  double gamma = gamma0;
  double sigma = sigma0;
  double beta = beta0;
  double t_sigma_gamma = t_sigma_gamma0;

  // Save objects initialization
  std::vector<MatIntCol> Xi_mcmc(niter, MatIntCol::Zero(H,Ttot));
  std::vector<MatCol> S_mcmc(niter, MatCol::Zero(H,Ttot));       
  std::vector<MatCol> U_mcmc(niter,  MatCol::Zero(H,Ttot));      
  std::vector<MatUnsCol> N_mcmc(niter,  MatUnsCol::Zero(Ttot,H));
  //std::vector<std::vector<MatUnsCol>> Dl_mcmc(niter, Dl0);
  std::vector<std::vector<MatCol>> Lambda_mcmc(niter, Lambda0);
  std::vector<std::vector<MatCol>> Lambda_star_mcmc(niter, std::vector<MatCol>(H));

  std::vector<double> phi_mcmc(niter,-1.0);   
  std::vector<double> gamma_mcmc(niter,-1.0); 
  std::vector<double> sigma_mcmc(niter,-1.0); 
  std::vector<double> beta_mcmc(niter,-1.0);  
  std::vector<double> t_sigma_gamma_mcmc(niter,-1.0);  

  // Usefull quantities for adaptive MCMC
  double s_adp = 2.83;
  MatCol Cadp = MatCol::Identity(2,2);
  MatCol Sigma_prop = var_beta * MatCol::Identity(2,2);
  VecCol RunningMean = VecCol::Zero(2);
  RunningMean(0) = gamma_mcmc[0]; RunningMean(1) = sigma_mcmc[0];
  MatCol gs = MatCol::Zero( niter_tot, 2 );
  double acc_prob = 0.0;

  // Start MCMC loop
  Rcpp::Rcout<<" Preprocessing finished. Start MCMC ... "<<std::endl;
  Progress progress_bar(niter_tot, print); // Initialize progress bar
  for(int it = 0; it < niter_tot; it++){

  	// ----------------------------------
  	//Rcpp::Rcout<<"UpdateXi"<<std::endl;
  	if(UpdateXi){
  		Xi = sample_Xi_tl(engine, Xi, S, N, phi, sigma, beta, t_sigma_gamma);
  	}
		// ----------------------------------
		//Rcpp::Rcout<<"UpdateS"<<std::endl;
  	if(UpdateS){
			S = sample_Stl(engine, Xi, U, phi, sigma, beta);
  	}
  	// ----------------------------------
  	//Rcpp::Rcout<<"UpdateU"<<std::endl;
  	if(UpdateU){
			U = sample_Utl(engine, S, t_sigma_gamma);
  	}
  	// ----------------------------------
  	//Rcpp::Rcout<<"UpdateLambda"<<std::endl;
  	if(UpdateLambda){
  		Lambda = sample_Lambda_itl(engine, Dl, Xi, delta);
  	}
  			// Check Lambda --> Is this necessary??
  			//for(int l=0; l < H; l++){
  				//VecCol aux = Lambda[l].colwise().sum();
  				//if (!((aux.array() - 1.0).abs() <= 1e-10).all()) {
  		    		//throw std::runtime_error("Error in Lambda: column values must sum to 1");
  				//}
  			//}
  	// ----------------------------------
  	//Rcpp::Rcout<<"UpdateDitl"<<std::endl;
  	if(UpdateDitl){
  		auto aux = sample_Ditl(engine, Lambda, Xi, D);
			Dl = aux.first;
			N  = aux.second;
  	}
  	// ----------------------------------
  	if(UpdatePhi || UpdateBeta || UpdateGamma || UpdateSigma){

		  	if(it > 100 ){
		  		Sigma_prop = s_adp*Cadp + 1e-6 * MatCol::Identity(2, 2);
		  	}

		  	//Rcpp::Rcout<<" +++++++++++++++++++++ "<<std::endl;
		  	//Rcpp::Rcout<<" it = "<<it<<std::endl;
		    //Rcpp::Rcout<<"Sigma_prop:"<<std::endl<<Sigma_prop<<std::endl;

  		VecCol aux = sample_hyparams_general( engine, Xi, S, phi, gamma, sigma, beta, t_sigma_gamma,
		                         								a_phi, b_phi, a_gamma, b_gamma, a_sigma, b_sigma, a_beta, b_beta,  
		                         								var_beta, Sigma_prop, s_adp, 
		                         								JointAdp, UpdatePhi, UpdateGamma, UpdateSigma, UpdateBeta);

  		//RunningMean(0) = ( (double)(it-1) )/( (double)(it) )*RunningMean(0) + 1.0/( (double)(it) ) * aux[1]; 
  		//RunningMean(1) = ( (double)(it-1) )/( (double)(it) )*RunningMean(1) + 1.0/( (double)(it) ) * aux[2]; 


	  	phi   = aux[0];
	  	gamma = aux[1]; 
	  	sigma = aux[2]; 
	  	beta  = aux[3]; 
	  	t_sigma_gamma = aux[4];
	  	acc_prob = aux[5];
  	}
  	// Update gs matrix and update adaptive hyperparameters
  	gs(it,0) = gamma; gs(it,1) = sigma;
  	if(it > 2 && it < 50000 && JointAdp){
			MatCol gs_it = gs.topRows(it);
			Cadp = my_cov( gs );
  		s_adp = std::exp( std::log(s_adp) + std::pow(it,-0.7)*( acc_prob - 0.234 ) );
  		if(s_adp > 10)
  			s_adp = 1.0;
  		if(s_adp < 1e-30)
  			s_adp = 1e-5;
  	}

  	// Save current iteration, if needed
  	if(it>=nburn && (it-nburn)%thin == 0){
  		if(nsaved >= niter)
  			throw std::runtime_error("Error, too many saved objects ");
  		 // Save objects
  		Xi_mcmc[nsaved] = Xi;
  		S_mcmc[nsaved]  = S;
  		U_mcmc[nsaved]  = U;
  		N_mcmc[nsaved]  = N;
  		//Dl_mcmc[nsaved] = Dl;
  		Lambda_mcmc[nsaved] = Lambda;

  		// compute and save Lambda_star
  		Rcpp::Rcout<<"Salvo Lambda_star_mcmc ... ";
  		for(int l=0; l < H; l++){
  		  std::vector<int> idx_born;  // vector with indexes when a trait is born
  		  std::vector<int> idx_surv;  // vector with indexes when a trait is survived
  		  std::vector<int> idx_noact; // vector with indexes when a trait is not active

  		  VecIntCol Xi_l = Xi.row(l);
  		  find_indices(Xi.row(l),idx_born,idx_surv,idx_noact);
  		  Lambda_star_mcmc[nsaved][l] = MatCol::Zero(idx_born.size(), V );
  		  for(int kk = 0; kk < idx_born.size(); kk++){
  		  	Lambda_star_mcmc[nsaved][l].row(kk) = Lambda[l].col(idx_born[kk]);	
  		  }
  		  
  		}
			Rcpp::Rcout<<" done! "<<std::endl;

			// Save hyperparameters
  		phi_mcmc[nsaved] = phi;
  		gamma_mcmc[nsaved] = gamma;
  		sigma_mcmc[nsaved] = sigma;
  		beta_mcmc[nsaved] = beta;
  		t_sigma_gamma_mcmc[nsaved] = t_sigma_gamma;

  		nsaved++;
    	//Rcpp::Rcout<<"it = "<<it<<std::endl;
    }

  	//throw std::runtime_error("FERMO IO ");
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

  Rcpp::Rcout<<"MCMC completed. Create return object ... "<<std::endl;
  return Rcpp::List::create(  
  	Rcpp::Named("Xi") = Xi_mcmc,
  	Rcpp::Named("S") = S_mcmc,
  	Rcpp::Named("U") = U_mcmc,
  	Rcpp::Named("N") = N_mcmc,
  	//Rcpp::Named("Dl") = Dl_mcmc,
  	Rcpp::Named("Lambda") = Lambda_mcmc,
  	Rcpp::Named("Lambda_star_mcmc") = Lambda_star_mcmc,
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


