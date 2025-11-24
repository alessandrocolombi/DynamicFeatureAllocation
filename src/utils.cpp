#include "utils.h"

// This function computes the inverse of a positive definite matrix.
// If the matrix is 2x2, the exact formula is used, hence in this case the matrix is not required to be pos. def
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
  my_posdef_inverse(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& A)
{
  // typedef
  using MatCol          = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using VecCol          = Eigen::VectorXd;

  if(A.rows() != A.cols())
    throw std::runtime_error("Error in my_posdef_inverse: matrix is not squared ");

  const unsigned int p{A.rows()};
  if(p == 1)
    return MatCol::Constant(1,1,1/A(0,0));
  else if(p == 2){
    MatCol res{MatCol::Constant(2,2,0.0)};
    res(0,0) = A(1,1);
    res(1,1) = A(0,0);
    res(0,1) = -A(0,1);
    res(1,0) = -A(1,0);
    return (res/(A(0,0)*A(1,1) - A(1,0)*A(0,1)));
  }
  else{
        Eigen::LLT<MatCol> lltOfA(A); // Compute the Cholesky decomposition
        // Check if the decomposition was successful
        if (lltOfA.info() == Eigen::NumericalIssue)
            throw std::runtime_error("Error in : Matrix is not positive definite. ");// decomposition failed
        MatCol I( MatCol::Identity(p,p) ); // define identity matrix
        return ( lltOfA.solve(I) ) ; // invert the matrix and return
  }
}


double my_var(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& data)
{
    int sz = data.size();
    if(sz == 0)
      throw std::runtime_error("Error in my_var: empty matrix inserted ");
    else if(sz == 1){
      return 0.0;
    }
    else{
      int n = data.rows();
      int p = data.cols();
      if( p > 1)
        throw std::runtime_error("Error in my_var: p is supposed to be 1 ");

      double var = (data.array() - data.mean()).array().square().sum()/(n-1.0);
      return var;
    }
}


Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
  my_cov(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& data)
{

    // typedef
    using MatCol          = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using VecCol          = Eigen::VectorXd;

    int sz = data.size();
    if(sz == 0)
      throw std::runtime_error("Error in my_cov: empty matrix inserted ");
    else{
      int n = data.rows();
      int p = data.cols();
      if( n == 1)
        return (MatCol::Constant(p,p,0.0));
      if( p == 1){
        double variance = my_var(data);
        return (MatCol::Constant(1,1,variance));
      }

      MatCol data_centered = data.rowwise() - data.colwise().mean(); // compute mean and center data matrix
      MatCol cov = (data_centered.transpose() * data_centered) / double(n - 1); // compute covariance matrix
      return cov;
    }

}

double log_dnorm(const double&x, const double& mean, const double& sd)
{
  double res{ -0.5*std::log(2.0*M_PI*sd*sd) - (0.5*(x-mean)*(x-mean) )/(sd*sd) };
  return res;
}

double log_dmvnorm(const Eigen::VectorXd& x,
                   const Eigen::VectorXd& mean,
                   const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& Sigma )
{
  const unsigned int p{x.size()};

  // checks
  if(p == 0)
    throw std::runtime_error("Error in log_dmvnorm: empty x inserted ");
  if(Sigma.rows() != Sigma.cols())
    throw std::runtime_error("Error in log_dmvnorm: Sigma is not squared ");
  if(Sigma.rows() != p )
    throw std::runtime_error("Error in log_dmvnorm: size of Sigma does not match with x size");
  if(mean.size() != p )
    throw std::runtime_error("Error in log_dmvnorm: size of mean does not match with x size");


  // 1-dim case
  if(p == 1)
    return log_dnorm(x[0], mean[0], Sigma(0,0));

  // p > 1 case:
  double res{ -((double)p/2.0)*std::log(2.0*M_PI) - 0.5*std::log(Sigma.determinant()) - 0.5 * (x-mean).transpose()* my_posdef_inverse(Sigma) *(x-mean) };
  return (res);
}


double log_dnct(const double&x, const double& dof, const double& location, const double& scale)
{
  // checks
  if(dof <= 0)
    throw std::runtime_error("Error in log_dnct: dof must be strictly positive ");
  if(scale <= 0)
    throw std::runtime_error("Error in log_dnct: scale must be strictly positive ");
  // Central student-t case is recovered setting location = 0.0 and scale = 1.0
  double res{0.0};
  res += gsl_sf_lngamma( (dof + 1.0)/2.0 ) - gsl_sf_lngamma( dof/2.0 ) - 0.5*std::log(M_PI*dof*scale*scale);
  res -= (dof + 1.0)/(2.0) * std::log( 1.0 + 1.0/dof * (x-location)/scale * (x-location)/scale );
  return (res);
}


double log_dmvt(const Eigen::VectorXd& x,
                const double& dof,
                const Eigen::VectorXd& location,
                const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& Scale )
{
  const unsigned int p{x.size()};

  // checks
  if(dof <= 0)
    throw std::runtime_error("Error in log_dnct: dof must be strictly positive ");
  if(p == 0)
    throw std::runtime_error("Error in log_dmvt: empty x inserted ");
  if(Scale.rows() != Scale.cols())
    throw std::runtime_error("Error in log_dmvt: Scale is not squared ");
  if(Scale.rows() != p )
    throw std::runtime_error("Error in log_dmvt: size of Scale does not match with x size");
  if(location.size() != p )
    throw std::runtime_error("Error in log_dmvt: size of location does not match with x size");


  // 1-dim case
  if(p == 1)
    return log_dnct(x[0], dof, location[0], Scale(0,0));

  // p > 1 case:
  double res{0.0};
  res += gsl_sf_lngamma( (dof + (double)p)/2.0 ) - gsl_sf_lngamma( dof/2.0 ) - ((double)p/2.0) * std::log(M_PI*dof) - 0.5*std::log(Scale.determinant());
  res -= (dof + (double)p)/(2.0) * std::log( 1.0 + 1.0/dof * (x-location).transpose()* my_posdef_inverse(Scale) *(x-location) );
  return (res);
}



std::vector<bool> find_all_int(const Rcpp::IntegerVector& vec,  const int& x)
{
  std::vector<bool> res;
  std::transform(  vec.cbegin(), vec.cend(), std::back_inserter(res), [& x](int vec_j) { return x == vec_j; }  );
  return res;
}

Rcpp::IntegerVector my_which_int(const Rcpp::IntegerVector& vec,  const int& x)
{
  int n = vec.size();
  std::vector<bool> idx{find_all_int(vec,x)};
  int nn = std::accumulate(idx.cbegin(), idx.cend(), 0);
  Rcpp::IntegerVector res(nn);
  if(nn == 0)
    return res;

  int ii = 0;
  for(std::size_t i = 0; i < n; i++){
    if(idx[i]){
      res[ii] = i;
      ii++;
    }
  }
  return res;
}

Eigen::MatrixXd myRowSlicing(const Eigen::MatrixXd& Mat, const std::vector<bool>& idx_rows )
{
  int n{Mat.rows()};
  int p{Mat.cols()};
  if(idx_rows.size() != n)
    throw std::runtime_error("Error in myRowSlicing: idx_rows must be a vector of size n (Mat.rows())");
  int nn = std::accumulate(idx_rows.cbegin(), idx_rows.cend(), 0);
  Eigen::MatrixXd res{Eigen::MatrixXd::Constant(nn,p,0.0)};
  if(nn == 0)
    return res;

  int ii = 0;
  for(std::size_t i = 0; i < n; i++){
    if(idx_rows[i]){
      res.row(ii) = Mat.row(i);
      ii++;
    }
  }
  return res;
}

std::map<int,int> table_Rcpp(const Rcpp::IntegerVector& s)
{
    std::map<int, int> res;
    int n = s.size();
    for(std::size_t i = 0; i < n; i++)
      res[ s[i] ]++;

    return res;
}



List classify_indices_cpp(const NumericVector& Z) 
{
  int Ttot = Z.size();

  std::vector<int> idx_born;
  std::vector<int> idx_surv;
  std::vector<int> idx_noact;

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

  return List::create(
    _["idx_born"] = idx_born,
    _["idx_surv"] = idx_surv,
    _["idx_noact"] = idx_noact
  );
}


double logZ_BFRY(double beta, double sigma, double t_const)
{
    const double EPS_TINY = 1e-16; // tolerance for tiny rounding errors
    
    if (beta <= 0)
        throw std::runtime_error("Error in logZ_BFRY: beta must be positive");
    if (t_const <= 0)
        throw std::runtime_error("Error in logZ_BFRY: t_const must be positive");

    if (sigma == 0 || sigma >= 1.0)
        throw std::runtime_error("Error in logZ_BFRY: sigma must be < 1 and != 0");

    // x = sigma * log(1 + t/beta)
    const double t_over_beta = t_const / beta;
    const double L = gsl_log1p(t_over_beta);   // safe for small t/beta
    const double x = sigma * L; 

    // Prepare safe expm1 value
    // For sigma > 0 we will need expm1(x) (>0).
    // For sigma < 0 we will need -expm1(x) (>0).
    double e = std::expm1(x); // exp(x) - 1
    double lg;
    if(sigma > 0.0){
      lg = std::lgamma(1.0 - sigma) - std::log(sigma); //lG(-sigma) = lG(sigma) - log(sigma)
    }
    else{
      e = -e;  // sigma < 0 : x < 0 so e = expm1(x) < 0 and -e = -expm1(x) > 0
      lg = std::lgamma(-sigma); //lG(-sigma)
    }
    // From this point on, e is supposed to be positive
    if(e <= 0.0) {
      // tolerate tiny negative rounding; fallback to linear approx
      if(e > -1e-5) {
        // when x ~ 0, expm1(x) ~ x
        e = (std::abs(x) > EPS_TINY ? std::abs(x) : EPS_TINY);
      } 
      else{
        Rcpp::Rcout<<"e = "<<e<<std::endl;
        Rcpp::Rcout<<"x = "<<x<<std::endl;
        Rcpp::Rcout<<"L = "<<L<<std::endl;
        Rcpp::Rcout<<"sigma = "<<sigma<<std::endl;
        Rcpp::Rcout<<"t_const = "<<t_const<<std::endl;
        Rcpp::Rcout<<"beta = "<<beta<<std::endl;
        throw std::runtime_error("Error in logZ_BFRY: unexpected nonpositive e");
      } 
    }
    double res = lg + sigma * std::log(beta) + std::log(e);    
    if(std::isnan(res))
       throw std::runtime_error("Error in logZ_BFRY: nan in res for sigma > 0");
    if (!std::isfinite(res)) {
       throw std::runtime_error("Error in logZ_BFRY: res is infinite");
    }
    return res; // return
}
