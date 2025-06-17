#ifndef __UTILS__H
#define __UTILS__H


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf.h>     //For special functions such as factorials

// Include file with basic libraries to include
#include "headers.h"


using namespace Rcpp;

// --------------------------------------------------------------------------------------------
// Rcpp implementation of some utility functions
// --------------------------------------------------------------------------------------------


// This function computes the inverse of a positive definite matrix.
// If the matrix is 2x2, the exact formula is used, hence in this case the matrix is not required to be pos. def
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
  my_posdef_inverse(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& A);


double my_var(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& data);

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
  my_cov(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& data);


// [[Rcpp::export]]
double log_dnorm(const double&x, const double& mean, const double& sd);

// [[Rcpp::export]]
double log_dmvnorm(const Eigen::VectorXd& x,
                   const Eigen::VectorXd& mean,
                   const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& Sigma );


// [[Rcpp::export]]
double log_dnct(const double&x, const double& dof, const double& location, const double& scale);

// [[Rcpp::export]]
double log_dmvt(const Eigen::VectorXd& x,
                const double& dof,
                const Eigen::VectorXd& location,
                const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& Scale );



std::vector<bool> find_all_int(const Rcpp::IntegerVector& vec,  const int& x);
Rcpp::IntegerVector my_which_int(const Rcpp::IntegerVector& vec,  const int& x);
Eigen::MatrixXd myRowSlicing(const Eigen::MatrixXd& Mat, const std::vector<bool>& idx_rows );
// [[Rcpp::export]]
std::map<int,int> table_Rcpp(const Rcpp::IntegerVector& s);


#endif
