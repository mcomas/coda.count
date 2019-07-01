// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat c_dirichlet_alr_approx(arma::vec alpha){
  int d = alpha.n_elem-1;
  //arma::vec mu(d);
  arma::mat pars(d,d+1);

  pars.head_cols(d).fill(R::trigamma(alpha[d]));  // SIGMA
  pars.col(d).fill(-R::digamma(alpha[d]));        // MU

  for(int i=0; i<d; i++){
    pars(i,i) += R::trigamma(alpha[i]);
    pars(i,d) += R::digamma(alpha[i]);
  }
  return pars;
}

// [[Rcpp::export]]
arma::mat c_logistic_ilr_approximation(int d){
  arma::mat res = arma::zeros(d, d+1);

  double ln_Mx = -(d+0.5) * log(d+1);
  for(int i=2; i < d+1; i++){
    ln_Mx += log(i);
  }
  res.diag().fill(exp(-2/d * ln_Mx - log(2*M_PI)));

  return(res);
}

// [[Rcpp::export]]
arma::mat c_logistic_alr_approximation(int d){
  arma::mat res(d, d+1);
  res.head_cols(d).fill(1);
  res.col(d).fill(0);
  res.diag() += 1;

  double ln_Mx = -(d+0.5) * log(d+1);
  for(int i=2; i < d+1; i++){
    ln_Mx += log(i);
  }
  res.head_cols(d) *= exp(-2/d * ln_Mx - log(2*M_PI));

  return(res);
}

// [[Rcpp::export]]
arma::mat c_gaussian_product(arma::mat pars1, arma::mat pars2){
  int d = pars1.n_rows;
  arma::mat res(d,d+1);

  arma::mat invSIGMA12 = arma::inv_sympd(pars1.head_cols(d) + pars2.head_cols(d));
  res.head_cols(d) = pars1.head_cols(d) * invSIGMA12 * pars2.head_cols(d);
  res.col(d) = pars2.head_cols(d) * invSIGMA12 * pars1.col(d) +
    pars1.head_cols(d) * invSIGMA12 * pars2.col(d);
  return(res);
}

// [[Rcpp::export]]
arma::mat c_gaussian_division(arma::mat pars1, arma::mat pars_res){
  int d = pars1.n_rows;
  arma::mat res(d,d+1);

  arma::mat invSIGMA_res = arma::inv_sympd(pars_res.head_cols(d));
  arma::mat invSIGMA1 = arma::inv_sympd(pars1.head_cols(d));

  res.head_cols(d) = arma::inv_sympd(invSIGMA_res - invSIGMA1);
  res.col(d) = res.head_cols(d) * (invSIGMA_res * pars_res.col(d) - invSIGMA1 * pars1.col(d));

  return(res);
}
