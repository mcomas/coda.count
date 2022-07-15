// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::mat pinv_sympd(arma::mat X){
  return arma::inv_sympd(X, arma::inv_opts::allow_approx);
}

// [[Rcpp::export]]
arma::mat pinv(arma::mat X){
  return arma::pinv(X);
}

int simplex_lattice_elements(int K, int SIZE){
  return R::choose(K+SIZE-1,SIZE);
}

// [[Rcpp::export]]
NumericMatrix c_simplex_lattice(int K, int SIZE) {
  int nrow = simplex_lattice_elements(K,SIZE);
  int k = K - 1;
  int k1 = k - 1;
  NumericMatrix out(nrow, K);
  NumericVector x(K);
  x(0) = SIZE;

  int target = 0;
  int i = -1;
  do{
    i++;
    out(i,_) = x;
    x(target) = x(target) - 1;
    if(target < k1){
      target = target + 1;
      x(target) = 1 + x(k);
      x(k) = 0;
    }else{
      x(k) = x(k) + 1;
      while(x(target) == 0){
        target = target - 1;
        if(target == -1){
          i++;
          out(i,_) = x;
          return(out);
        }
      }
    }
  }while(true);
  return(out);
}

