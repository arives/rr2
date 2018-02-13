// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' A loop in c++ to speed up
//'
//' @param R a matrix (not sparse matrix, use \code{as.matrix(R)} to convert if necessary)
//' @param V a symmetric matrix (not sparse matrix)
//' @export
// [[Rcpp::export]]
arma::vec loop_cpp(const arma::vec& R, const arma::mat& V) {
  // Initiate Rhat
  int n = R.size();
  arma::vec Rhat(n); // to save results
  arma::uvec idx(n); // to identify elements to remove
  idx.fill(1.0);
  
  for (int j = 0; j < n; j++) {
    idx(j) = 0.0;
    arma::uvec idx2 = find(idx);
    arma::vec r = R.elem(idx2); // R[-j]
    arma::mat VV = V.submat(idx2, idx2); // V[-j, -j]
    arma::mat iVV = arma::inv_sympd(VV);
    rowvec v1 = V.row(j);  // V[j, ]
    vec v = v1.elem(idx2); // V[j, -j]
    double res = as_scalar(trans(v) * iVV * r);
    Rhat(j) = res;
    idx.fill(1);
  }
  return Rhat;
}

/***R
# loop_cpp(as.matrix(R), as.matrix(V))
*/
