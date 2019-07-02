#include "coda_base.h"

arma::mat inv_B_coordinates(arma::mat ilrX, arma::mat Binv){

  arma::mat X = arma::exp(ilrX * Binv);
  // Closuring
  arma::mat sumX = arma::sum(X, 1);
  for(unsigned int i = 0; i < X.n_cols; i++){
    X.col(i) = X.col(i) / sumX;
  }
  return X;
}

// [[Rcpp::export]]
arma::mat H(int d){
  return(arma::eye(d,d) + arma::ones(d,d));
}

// [[Rcpp::export]]
arma::mat F(int d){
  arma::mat res = arma::zeros(d,d+1);
  res.diag().fill(1);
  res.col(d).fill(-1);
  return(res);
}
