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

// [[Rcpp::export]]
arma::mat ilr_basis_default(unsigned int dim){
  arma::mat B = arma::mat(dim, dim-1);

  for(unsigned int i = 0; i < dim - 1; i++){
    unsigned int I1 = i + 1;
    unsigned int I2 = i + 2;
    double l = 1/std::sqrt((double)(I1*I2));
    double r = - std::sqrt((double)I1/I2);
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
