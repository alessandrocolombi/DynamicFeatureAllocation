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


int sample_Ditl(sample::GSL_RNG const & engine, const std::vector<MatCol>& Lambda_itl, const MatIntCol& Xi, const MatIntCol& D) 
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

  return 1;
}


// It should be std::vector<MatUnsCol>& D_itl
// but for testing purposes I am using std::vector<MatCol>& D_itl
int sample_Lambda_itl(sample::GSL_RNG const & engine, const std::vector<MatCol>& D_itl, const MatIntCol& Xi, const double& delta)
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

  return 1;
}

int sample_Stl(sample::GSL_RNG const & engine, const MatIntCol& Xi, const MatCol& U, const double& phi, const double& sigma, const double& b)
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
      double shape = Xi(l,t) - sigma;
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
  return 1;
}

int sample_Utl(sample::GSL_RNG const & engine, const MatCol& S, const double& t_sigma_gamma)
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
      double Y = runif(engine);
      double temp = invCDF(S(l,t), Y);
      if( temp <= 0 || std::isnan(temp) ){
        Rcpp::Rcout<<"t = "<<t<<"; l = "<<l<<"temp = "<<temp<<std::endl;
        throw std::runtime_error("Error in sample_Utl: U is out of range ");
      }
    }
  }

  return 1;

}

int sample_Xi_tl(sample::GSL_RNG const & engine, const MatIntCol& Xi_old, const MatCol& S, const MatUnsCol& N_tl, 
                 const double& phi, const double& sigma, const double& b, const double& t_sigma_gamma)
{
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

  for(int l=0; l < H; l++){
    for(int t=0; t < Ttot; t++){
      int xi_prime = rpoisson(phi * S(l,t));
      int xi_tl = Xi_old(l,t);
      double log_R_tl = (double)xi_tl - (double)xi_prime + N_tl(t,l) * ( std::log((double)xi_prime) - std::log((double)xi_tl) );
      if(t < Ttot - 1){
        log_R_tl += std::log((double)xi_prime - sigma) - std::log((double)xi_tl - sigma) + std::lgamma(1.0 + (double)xi_tl - sigma) - std::lgamma(1.0 + (double)xi_prime - sigma);
        log_R_tl += ((double)xi_prime - (double)xi_tl)*std::log( S(l,t+1) );
        log_R_tl += std::log( std::pow(b+phi+t_sigma_gamma,sigma-(double)xi_tl) - std::pow(b+phi,sigma-(double)xi_tl) ) - std::log( std::pow(b+phi+t_sigma_gamma,sigma-(double)xi_prime) - std::pow(b+phi,sigma-(double)xi_prime) );
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

    }
  }
  return 1;
}