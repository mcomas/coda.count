#include "coda_base.h"


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
