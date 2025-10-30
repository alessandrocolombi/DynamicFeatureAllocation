#include "FC.h"

using namespace Rcpp;


void find_indices(const VecIntCol& Z, std::vector<int>& idx_born, std::vector<int>& idx_surv, std::vector<int>& idx_noact)
{
  int Ttot = Z.size();
  // loop through elements
  for (int t = 0; t < Ttot; t++) {
    if (Z[t] > 0) {
      // If here, the feature at time t is active
      if (t == 0 || Z[t - 1] == 0) {
        // born: first element = 1, or preceded by 0
        idx_born.push_back(t);
      } else {
        // 1 not preceded by 0 → survivor
        idx_surv.push_back(t);
      }
    } else {
      // zero → no activity
      idx_noact.push_back(t);
    }
  }
}


std::pair<std::vector<MatUnsCol>, MatUnsCol> 
sample_Ditl(sample::GSL_RNG const & engine, const std::vector<MatCol>& Lambda_itl, const MatIntCol& Xi, const MatIntCol& D) 
{
  sample::rmultinomial<VecUnsCol> rmultinomial; // define callable object to generate random samples from a Multivariate distribution

  const int H = Xi.rows(); // Number of atoms
  const int Ttot = Xi.cols(); // Time windomw
  if(Lambda_itl.size() != H)
    throw std::runtime_error("Error in sample_Ditl: Xi.rows() != Lambda_itl.size()");
  if(H <= 0)
    throw std::runtime_error("Error in sample_Ditl: H must be >= 1 ");
  if(Lambda_itl[0].cols() != Ttot)
    throw std::runtime_error("Error in sample_Ditl: Xi.cols() != Lambda_itl[0].cols()");

  const int V = Lambda_itl[0].rows(); // Vocabulary size
  
  
  std::vector<MatUnsCol> D_itl(H, MatUnsCol::Zero(V, Ttot)); // Define main object
    
  for(int i = 0; i < V; i++) {
    for(int t = 0; t < Ttot; t++) {
      VecCol zeta_it{VecCol::Zero(H)}; // inizialize weights to 0
      // Build zeta_it = Lambda[l](i,t) * Xi(l,t)
      for(int l = 0; l < H; l++) {
        zeta_it(l) = Lambda_itl[l](i,t) * Xi(l,t);
        if(std::isnan(zeta_it(l)))
          throw std::runtime_error("Error in sample_Ditl, get a nan in zeta_it ");
      }
      double sum_w = zeta_it.sum(); // compute the sum
      if(sum_w > 0) { // if 0, do nothing. otherwise ..
        zeta_it /= sum_w;  // normalize
        //Rcpp::Rcout<<"zeta_it:"<<std::endl<<zeta_it<<std::endl;
        int n = D(i,t); // number of trials
        VecUnsCol temp = rmultinomial(engine, n, zeta_it); // sample from multinomial distribution
        for(int l = 0; l < H; l++) { 
          D_itl[l](i,t) = temp(l); // Save values in D_itl
        }
      
        //Rcpp::Rcout<<"temp: "<<temp<<std::endl;
      }
      // End of time t for word i
    }
    // End of word i for all t
  }
  // End of all words, all times
  MatUnsCol N_tl{MatUnsCol::Zero(Ttot,H)}; // Define matrix with total counts for each topic
  for(int l = 0; l < H; l++) { 
    //Rcpp::Rcout<<"D_itl["<<l<<"]:"<<std::endl<<D_itl[l]<<std::endl;
    N_tl.col(l) = D_itl[l].colwise().sum();
  }
  //Rcpp::Rcout<<"N_tl:"<<std::endl<<N_tl<<std::endl;

  return std::make_pair(D_itl, N_tl);
}

std::vector<MatCol> sample_Lambda_itl(sample::GSL_RNG const & engine, const std::vector<MatUnsCol>& D_itl, const MatIntCol& Xi, const double& delta)
{
  sample::rgamma rgamma; // define callable object to generate random samples from a gamma distribution

  const int H = Xi.rows(); // Number of atoms
  const int Ttot = Xi.cols(); // Time windomw
  if(D_itl.size() != H)
    throw std::runtime_error("Error in sample_Ditl: Xi.rows() != D_itl.size()");
  if(H <= 0)
    throw std::runtime_error("Error in sample_Ditl: H must be >= 1 ");
  if(D_itl[0].cols() != Ttot)
    throw std::runtime_error("Error in sample_Ditl: Xi.cols() != D_itl[0].cols()");
  if(delta <= 0)
    throw std::runtime_error("Error in sample_Ditl: delts must be positive ");

  const int V = D_itl[0].rows(); // Vocabulary size

  std::vector<MatCol> Lambda_itl(H, MatCol::Zero(V, Ttot)); // Define main object

  for(int l=0; l < H; l++){
      std::vector<int> idx_born;  // vector with indexes when a trait is born
      std::vector<int> idx_surv;  // vector with indexes when a trait is survived
      std::vector<int> idx_noact; // vector with indexes when a trait is not active

      VecIntCol Xi_l = Xi.row(l);
      find_indices(Xi.row(l),idx_born,idx_surv,idx_noact);

      // These features must be updated from a Dirichlet with updated hyperparameters
      for(int t : idx_born){
        VecCol temp{VecCol::Zero(V)}; 
        for(int i = 0; i < V; i++){
          double shape = delta;
          int s = t; // start from time t
          bool flag = TRUE; // go on with the sum
          while(flag & s < Ttot){
            shape += (double)D_itl[l](i,s);
            if(s+1 >= Ttot || Xi_l[s+1] == 0)
              flag = FALSE; // if so, stop
            s++;
          }
          temp(i) = rgamma(engine, shape, 1.0);
        }
        temp /= temp.sum(); // normalize to get a Dirichlet distribution
        Lambda_itl[l].col(t) = temp ; // save new value from posterior
      }
      // These features are deterministically assigned
      for (int t : idx_surv){
        if(t == 0)
          throw std::runtime_error("Error in sample_Lambda_itl: Features are time t=0 can not be survival features");
        Lambda_itl[l].col(t) = Lambda_itl[l].col(t-1);
      }
      for (int t : idx_noact){
        VecCol temp{VecCol::Zero(V)}; 
        for(int i = 0; i < V; i++){
          temp(i) = rgamma(engine, delta, 1.0); // no update from the data
        }
        temp /= temp.sum(); // normalize to get a Dirichlet distribution
        Lambda_itl[l].col(t) = temp ; // save draw from the prior
      }
  }

  return Lambda_itl;
}

MatCol sample_Stl(sample::GSL_RNG const & engine, const MatIntCol& Xi, const MatCol& U, const double& phi, const double& sigma, const double& b)
{
  sample::rgamma rgamma; // define callable object to generate random samples from a gamma distribution

  const int H = Xi.rows(); // Number of atoms
  const int Ttot = Xi.cols(); // Time window
  if(H <= 0)
    throw std::runtime_error("Error in sample_Stl: H must be >= 1 ");
  if(U.rows() != H)
    throw std::runtime_error("Error in sample_Stl: U.rows() != H ");
  if(U.cols() != Ttot)
    throw std::runtime_error("Error in sample_Stl: U.cols() != Ttot ");
  if( phi <= 0 || b <= 0 )
    throw std::runtime_error("Error in sample_Stl: phi or b are null or negative ");
  if( sigma < 0 || sigma >= 1)
    throw std::runtime_error("Error in sample_Stl: sigma is out of range");

  MatCol S{ MatCol::Zero(H,Ttot) }; // inizialize main object
  for(int l=0; l < H; l++){
    for(int t=0; t < Ttot; t++){
      double shape = 1.0 + Xi(l,t) - sigma;
      double rate = U(l,t) + b + phi;
      if(t > 0){
        shape += Xi(l,t-1);
        rate += phi;
      }
      if( shape <= 0 || rate <= 0 || std::isnan(shape) || std::isnan(rate) )
          throw std::runtime_error("Error in sample_Stl, invalid shape or rate");
      S(l,t) = rgamma(engine, shape, 1.0/rate);
    }
  }
  return S;
}

MatCol sample_Utl(sample::GSL_RNG const & engine, const MatCol& S, const double& t_sigma_gamma)
{
  sample::runif runif; // define callable object to generate random samples from a Uniform distribution

  const int H = S.rows(); // Number of atoms
  const int Ttot = S.cols(); // Time window
  if(H <= 0)
    throw std::runtime_error("Error in sample_Utl: H must be >= 1 ");
  if( t_sigma_gamma <= 0 || std::isnan(t_sigma_gamma) )
    throw std::runtime_error("Error in sample_Utl: t_sigma_gamma is out of range ");

  MatCol U{ MatCol::Zero(H,Ttot) }; // inizialize main object

  auto invCDF = [t_sigma_gamma](double s, double y){
    double temp = - (y * ( 1.0 - exp(-s*t_sigma_gamma) ));
    temp = gsl_log1p( temp );
    temp = - (1.0/s)*temp;
    return temp;
  };

  for(int l=0; l < H; l++){
    for(int t=0; t < Ttot; t++){
      if(S(l,t) <= 0)
        throw std::runtime_error("Error in sample_Utl: S_tl can not be negative or zero ");
      
      double Y = runif(engine);
      double temp = invCDF(S(l,t), Y);
      if( temp <= 0 || std::isnan(temp) ){
        Rcpp::Rcout<<"t = "<<t<<"; l = "<<l<<"temp = "<<temp<<std::endl;
        throw std::runtime_error("Error in sample_Utl: U is out of range ");
      }
      U(l,t) = temp;
    }
  }

  return U;
}

MatIntCol sample_Xi_tl(sample::GSL_RNG const & engine, const MatIntCol& Xi_old, const MatCol& S, const MatUnsCol& N_tl, 
                       const double& phi, const double& sigma, const double& b, const double& t_sigma_gamma)
{
  //Rcpp::Rcout<<"Dentro!"<<std::endl;

  sample::rpoisson rpoisson; // define callable object to generate random samples from a gamma distribution
  sample::runif runif; // define callable object to generate random samples from a Uniform distribution

  const int H = S.rows(); // Number of atoms
  const int Ttot = S.cols(); // Time window
  if(H <= 0)
    throw std::runtime_error("Error in sample_Xi_tl: H must be >= 1 ");
  if(Xi_old.cols() != Ttot)
    throw std::runtime_error("Error in sample_Xi_tl: Xi_old.cols() != Ttot ");
  if(Xi_old.rows() != H)
    throw std::runtime_error("Error in sample_Xi_tl: Xi_old.rows() != H ");
  if(N_tl.rows() != Ttot)
    throw std::runtime_error("Error in sample_Xi_tl: N_tl.rows() != Ttot ");
  if(N_tl.cols() != H)
    throw std::runtime_error("Error in sample_Xi_tl: N_tl.cols() != H ");
  if( phi <= 0 || b <= 0 )
    throw std::runtime_error("Error in sample_Xi_tl: phi or b are null or negative ");
  if( sigma < 0 || sigma >= 1)
    throw std::runtime_error("Error in sample_Xi_tl: sigma is out of range");
  if( t_sigma_gamma <= 0 || std::isnan(t_sigma_gamma) )
    throw std::runtime_error("Error in sample_Xi_tl: t_sigma_gamma is out of range ");

  MatIntCol Xi{MatIntCol::Zero(H,Ttot)}; // initialize main object
  double b_phi   = b+phi;
  double b_phi_t = b+phi+t_sigma_gamma;
  double temp    = std::log(b_phi) - std::log(b_phi_t);

  for(int l=0; l < H; l++){
    for(int t=0; t < Ttot; t++){
      if(S(l,t) <= 0)
        throw std::runtime_error("Error in sample_Xi_tl: S_tl can not be negative or zero ");

      int xi_prime = rpoisson(phi * S(l,t));
      int xi_tl = Xi_old(l,t);
      double log_R_tl =  + N_tl(t,l) * ( std::log((double)xi_prime) - std::log((double)xi_tl) );
      Rcpp::Rcout<<"log_R_tl 1 = "<<log_R_tl<<std::endl;
      if(t < Ttot - 1){
        log_R_tl += std::log((double)xi_prime - sigma) - std::log((double)xi_tl - sigma) + std::lgamma(1.0 + (double)xi_tl - sigma) - std::lgamma(1.0 + (double)xi_prime - sigma);
        Rcpp::Rcout<<"log_R_tl 2 = "<<log_R_tl<<std::endl;
        log_R_tl += ((double)xi_prime - (double)xi_tl)*std::log( S(l,t+1) );
        Rcpp::Rcout<<"log_R_tl 3 = "<<log_R_tl<<std::endl;
        log_R_tl += ( (double)xi_prime - (double)xi_tl ) * std::log(b_phi_t);
        Rcpp::Rcout<<"log_R_tl 4 = "<<log_R_tl<<std::endl;
        double x_old   = sigma - (double)xi_tl;
        double x_prime = sigma - (double)xi_prime;

        log_R_tl += std::log( -gsl_expm1( x_old  *temp ) );
        Rcpp::Rcout<<"log_R_tl 5 = "<<log_R_tl<<std::endl;
        log_R_tl -= std::log( -gsl_expm1( x_prime*temp ) );
        Rcpp::Rcout<<"log_R_tl 6 = "<<log_R_tl<<std::endl;

        //log_R_tl += std::log( std::pow(b+phi+t_sigma_gamma,sigma-(double)xi_tl) - std::pow(b+phi,sigma-(double)xi_tl) ) - std::log( std::pow(b+phi+t_sigma_gamma,sigma-(double)xi_prime) - std::pow(b+phi,sigma-(double)xi_prime) );
        
      }
      if( std::isnan(log_R_tl) ){
        Rcpp::Rcout<<"("<<l<<","<<t<<") : "<<xi_tl<<" vs "<<xi_prime<<", S(l,t) = "<<S(l,t)<<std::endl;
        throw std::runtime_error("Error in sample_Xi_tl: get nan in acceptance probability");
      }

      double u = runif(engine);
      if(std::log(u) < log_R_tl ){
        // accepted
        Xi(l,t) = xi_prime;
      }
      else{
        // rejected
        Xi(l,t) = xi_tl;
      }

      Rcpp::Rcout<<"("<<l<<","<<t<<") : "<<xi_tl<<" vs "<<xi_prime<<", prob "<<std::exp(log_R_tl)<<std::endl;
    }
    if(l == 4)
      throw std::runtime_error("FERMO IO ");
  }
  return Xi;
}


// This function implements the proposed update in Section B.4 of [C.Naik, F.Caron, J.Rousseau, Y.Teh, K.Palla] Bayesian Nonparametrics for Sparse Dynamic Networks (2023)
VecCol sample_hyparams( sample::GSL_RNG const & engine, const MatIntCol& Xi, const MatCol& S, 
                        const double& phi_old, const double& gamma_old, const double& sigma_old, const double& b_old, 
                        const double& t_sigma_gamma_old,
                        const double& a_phi, const double& b_phi, 
                        const double& a_gamma, const double& b_gamma, 
                        const double& a_sigma, const double& b_sigma, 
                        const double& a_beta, const double& b_beta, 
                        const double& var_phi, const double& var_gamma, const double& var_sigma, const double& var_beta,
                        bool UpdatePhi, bool UpdateGamma, bool UpdateSigma, bool UpdateBeta)
{
  sample::rnorm rnorm; // define callable object to generate random samples from a normal distribution
  sample::runif runif; // define callable object to generate random samples from a uniform distribution

  const int H = S.rows(); // Number of atoms
  const int Ttot = S.cols(); // Time window
  if(H <= 0)
    throw std::runtime_error("Error in sample_hyparams: H must be >= 1 ");
  if(Xi.cols() != Ttot)
    throw std::runtime_error("Error in sample_hyparams: Xi.cols() != Ttot ");
  if(Xi.rows() != H)
    throw std::runtime_error("Error in sample_hyparams: Xi.rows() != H ");
  if( phi_old <= 0 || b_old <= 0 || gamma_old <= 0)
    throw std::runtime_error("Error in sample_hyparams: phi_old, b_old or gamma_old are null or negative ");
  if( sigma_old < 0 || sigma_old >= 1)
    throw std::runtime_error("Error in sample_hyparams: sigma_old is out of range");
  if( t_sigma_gamma_old <= 0 || std::isnan(t_sigma_gamma_old) )
    throw std::runtime_error("Error in sample_hyparams: t_sigma_gamma is out of range ");

  // Sample proposed values
  double phi_prime{phi_old};
  double gamma_prime{gamma_old};
  double sigma_prime{sigma_old};
  double b_prime{b_old};
  double t_sigma_gamma_prime{t_sigma_gamma_old};
  
  double lognew; // auxiliary
  double logitnew; // auxiliary
  if(UpdatePhi){
    // Log-normal proposal
    lognew = rnorm( engine, std::log(phi_old), std::sqrt(var_phi) );
    phi_prime = std::exp(lognew);
  }
  if(UpdateGamma){
    // Log-normal proposal
    lognew = rnorm( engine, std::log(gamma_old), std::sqrt(var_gamma) );
    gamma_prime = std::exp(lognew);
  }
  if(UpdateSigma){
    // Logit-normal proposal
    double logit_old = std::log( sigma_old ) - std::log( 1.0 - sigma_old );
    logitnew = rnorm( engine, logit_old, std::sqrt(var_sigma) );
    sigma_prime = 1.0/( 1 + std::exp(-logitnew) );
  }
  if(UpdateBeta){
    // Log-normal proposal
    lognew = rnorm( engine, std::log(b_old), std::sqrt(var_beta) );
    phi_prime = std::exp(lognew);
  }
  if(UpdateGamma || UpdateSigma){
    // Update t_gamma_sigma
    t_sigma_gamma_prime = std::exp( 1.0/sigma_prime * ( std::log(sigma_prime*(double)H) - std::log(gamma_prime) )  );
  }

  // Compute acceptance prob.
  double log_prop_phi   = a_phi*  ( std::log(phi_prime)  -std::log(phi_old)    ) - b_phi*(phi_prime-phi_old);
  double log_prop_gamma = a_gamma*( std::log(gamma_prime)-std::log(gamma_old)  ) - b_gamma*(gamma_prime-gamma_old);
  double log_prop_sigma = a_sigma*( std::log(sigma_prime)-std::log(sigma_prime)) + b_sigma*(std::log(1.0 - sigma_prime)-std::log(1.0 - sigma_prime)) ;
  double log_prop_beta  = a_beta* ( std::log(b_prime)    -std::log(b_old)      ) - b_beta*(b_prime-b_old);

  // Auxiliary quantities
  double sumS0 = S.col(0).sum();
  double sumXi0 = Xi.col(0).sum();
  double sumS  = S.sum();
  double sumXi = Xi.sum();

  // Target ratio
  double log_target = (double)H * (std::log(sigma_prime) - std::log(sigma_old) + std::lgamma(1.0-sigma_old) - std::lgamma(1.0-sigma_prime) );
  log_target += (b_old+phi_old-b_prime-phi_prime)*sumS0;

  log_target += (double)H * ( sigma_old*std::log(b_old+t_sigma_gamma_old) + 
                              std::log( -gsl_expm1(sigma_old*( std::log(b_old) - std::log(b_old+t_sigma_gamma_old) ))  ) );
  log_target -= (double)H * ( sigma_prime*std::log(b_prime+t_sigma_gamma_prime) + 
                              std::log( -gsl_expm1(sigma_prime*( std::log(b_prime) - std::log(b_prime+t_sigma_gamma_prime) ))  ) );
  
  log_target += (sigma_old-sigma_prime)*sumS + (std::log(phi_prime)-std::log(phi_old))*sumXi;
  log_target += (b_old + 2*phi_old - b_prime - 2.0*phi_prime)*(sumS - sumS0);

  log_target += ( sigma_old  *(double)H*(double)(Ttot-1) - (sumXi-sumXi0) ) * std::log(b_old+phi_old+t_sigma_gamma_old);
  log_target -= ( sigma_prime*(double)H*(double)(Ttot-1) - (sumXi-sumXi0) ) * std::log(b_prime+phi_prime+t_sigma_gamma_prime);

  for(int t=0; t < Ttot; t++){
    for(int l=0; l < H; l++){
      double ynew = t_sigma_gamma_prime*S(l,t);
      double yold = t_sigma_gamma_old*S(l,t);
      log_target += std::log( -gsl_expm1(-ynew) ) - std::log( -gsl_expm1(-yold) );
      if(t > 0){
        log_target += std::lgamma( 1.0 - sigma_old + Xi(l,t-1) ) - std::lgamma( 1.0 - sigma_prime + Xi(l,t-1) );
        log_target += std::log( sigma_prime - Xi(l,t-1) ) - std::log( sigma_old - Xi(l,t-1) );
        // Final term in Equation (37) Naik:
        double A_old = b_old + phi_old;
        double t_old = t_sigma_gamma_old;
        double x_old = sigma_old - Xi(l,t-1);
        double A_prime = b_prime + phi_prime;
        double t_prime = t_sigma_gamma_prime;
        double x_prime = sigma_prime - Xi(l,t-1);

        log_target += std::log( -gsl_expm1( x_old*(std::log(A_old)     - std::log(A_old+t_old)) ) );
        log_target -= std::log( -gsl_expm1( x_prime*(std::log(A_prime) - std::log(A_prime+t_prime)) ) );
      }

    }
  }

  // Final log acceptance probability
  double log_acc = log_target + log_prop_phi + log_prop_gamma + log_prop_sigma + log_prop_beta;

  // Accept / Reject the move
  double phi_res{phi_old};
  double gamma_res{gamma_old};
  double sigma_res{sigma_old};
  double b_res{b_old};
  double t_sigma_gamma_res{t_sigma_gamma_old};

  double u = runif(engine);
  if( std::log(u) < log_acc ){
    phi_res = phi_prime;
    gamma_res = gamma_res;
    sigma_res = sigma_prime;
    b_res = b_prime;
    t_sigma_gamma_res = t_sigma_gamma_prime;
  }
  // Checks:
  if( phi_res <= 0 )
    throw std::runtime_error("Error in sample_hyparams: phi can not be negative or zero ");
  if( gamma_res <= 0 )
    throw std::runtime_error("Error in sample_hyparams: gamma can not be negative or zero ");
  if( b_res <= 0 )
    throw std::runtime_error("Error in sample_hyparams: b can not be negative or zero ");
  if( sigma_res <= 0.0 || sigma_res >= 1.0 )
    throw std::runtime_error("Error in sample_hyparams: sigma must be in (0,1) ");
  if( t_sigma_gamma_res <= 0.0 )
    throw std::runtime_error("Error in sample_hyparams: t_sigma_gamma_res can not be negative or zero ");
  
  // Return
  VecCol res_vec(5);
  res_vec << phi_res, gamma_res, sigma_res, b_res, t_sigma_gamma_res;
  return res_vec;
}