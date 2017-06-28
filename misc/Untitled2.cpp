#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat invcpp(arma::mat V){
  return arma::pinv(V);
}


// [[Rcpp::export]]
arma::vec rm1vec(arma::vec x, arma::uword j){
  int n = x.size();
  arma::uvec idx(n);
  idx.fill(1);
  idx(j - 1) = 0;
  return x.elem(find(idx));
}

// [[Rcpp::export]]
double rep1(arma::mat x, vec y){
  return as_scalar(x * y);
  // return x - y;
}

// [[Rcpp::export]]
vec subjj(arma::mat V, uword j){
  rowvec v = V.row(j - 1);
  Rcout << v << endl;
  vec idx(v.size());
  idx.fill(1);
  idx(j - 1) = 0;
  Rcout << idx << endl;
  return v.elem(find(idx));
  // return x - y;
}

// [[Rcpp::export]]
double onetime(arma::vec R, arma::mat V, arma::mat iV, arma::mat X) {
  // Initiate Rhat
  int n = R.size();
  arma::vec Rhat(n); // to save results
  arma::uvec idx(n);
  idx.fill(1);

  int j = 0;
    idx(j) = 0;
    arma::vec r = R.elem(find(idx)); // r[-j]
    arma::mat VV = V.submat(find(idx), find(idx)); // V[-j, -j]
    // arma::mat VV = sub1(V, j); // V[-j, -j]

    // arma::mat VV = V.shed_row(j);

    arma::mat iVV = arma::pinv(VV);

    arma::mat AA = X.t() * iV * X;
    arma::mat BB = X.t() * iV * R;
    arma::mat bbhat = arma::solve(AA, BB);
    bbhat.print();
    // arma::uvec first_elem(2);
    // first_elem.fill(0);
    //
    // arma::vec bbhat1 = bbhat.elem(first_elem);

    double bbhat1 = as_scalar(bbhat.row(0));
    Rcout << bbhat1 << endl;

    // arma::vec v = V.row(j);
    rowvec v1 = V.row(j);
    v1.print();
    mat v = v1.elem(find(idx));
    v.print();

    // arma::vec res = bbhat + v * iVV * (r - bbhat);
    double rsi = as_scalar(trans(v) * iVV * (r - bbhat1));
    Rcout << rsi << endl;
    double res = bbhat1 + rsi;

  return res;
}

/*** R
# as.matrix(solve(VV))
# invcpp(as.matrix(VV))
rm1vec(1:10, 2)
rep1(matrix(1:10, nrow = 1), 1:10)
subjj(matrix(1:16, 4, 4), 2)
onetime(as.matrix(R), as.matrix(V), as.matrix(iV), X)
*/
