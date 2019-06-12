// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "coda.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::mat ilr_basis(unsigned int dim){
  arma::mat B = arma::mat(dim, dim-1);

  for(unsigned int i = 0; i < dim - 1; i++){
    unsigned int I1 = i + 1;
    unsigned int I2 = i + 2;
    double l = 1/sqrt(I1*I2);
    double r = - sqrt((double)I1/I2);
    for(unsigned int j = 0; j < I1; j++){
      B(j,i) = l;
    }
    B(I1,i) = r;
    for(unsigned int j = I2; j < dim; j++){
      B(j,i) = 0;
    }
  }
  return(B);
}

//' @export
// [[Rcpp::export]]
arma::mat ilr_coordinates(arma::mat X){
  arma::mat B = ilr_basis(X.n_cols);

  arma::mat logX= log(X);

  return(logX * B);
}

//' @export
// [[Rcpp::export]]
arma::mat inv_ilr_coordinates(arma::mat ilrX){
  arma::mat B = ilr_basis(ilrX.n_cols + 1);

  arma::mat X = arma::exp(ilrX * B.t());
  // Closuring
  arma::mat sumX = arma::sum(X, 1);
  for(unsigned int i = 0; i < X.n_cols; i++){
    X.col(i) = X.col(i) / sumX;
  }
  return X;
}

//' @export
// [[Rcpp::export]]
arma::mat ilr_to_alr(arma::mat ILR){
  int k = ILR.n_cols;
  arma::mat B = ilr_basis(k+1);
  arma::mat ILR_TO_ALR = B(arma::span(0,k-1),arma::span(0,k-1));
  for(int i =0;i < k;i++){
    ILR_TO_ALR.col(i) = ILR_TO_ALR.col(i) - B(k,i);
  }
  return(ILR * ILR_TO_ALR.t());
}

//' @export
// [[Rcpp::export]]
arma::mat alr_to_ilr(arma::mat ALR){
  int k = ALR.n_cols;
  arma::mat B = ilr_basis(k+1);
  arma::mat ALR_TO_ILR = B(arma::span(0,k-1),arma::span(0,k-1));
  for(int i =0;i < k;i++){
    ALR_TO_ILR.col(i) = ALR_TO_ILR.col(i) - B(k,i);
  }
  return(ALR * ALR_TO_ILR.t().i());
}
