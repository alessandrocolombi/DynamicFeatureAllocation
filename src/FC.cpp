#include "FC.h"

using namespace Rcpp;


void find_indices(const VecIntCol& Z, std::vector<int>& idx_born, std::vector<int>& idx_surv, std::vector<int>& idx_noact)
{
  int Ttot = Z.size();
  // loop through elements
  for (int t = 0; t < Ttot; t++) {
    if (Z[t] > 0) {
      // If here, the feature at time t is active
      if (t == 0 || Z[t - 1] > 0) {
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


int sample_Ditl(sample::GSL_RNG const & engine,
				        const std::vector<MatCol>& Lambda_itl, 
				        const MatIntCol& Xi, 
				        const MatIntCol& D) 
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


int sample_Lambda_itl(sample::GSL_RNG const & engine, const std::vector<MatUnsCol>& D_itl, const MatIntCol& Xi, const double& delta)
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
}

