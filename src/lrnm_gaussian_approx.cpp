// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "gaussian_approx.h"
#include "lrnm_utils.h"
#include "coda_base.h"

//' @export
// [[Rcpp::export]]
double c_d_lrnm_gaussian_approx(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv){

  double integral = 0;

  return integral;
}

//' @export
// [[Rcpp::export]]
arma::mat c_posterior_alr_approximation(arma::vec x, arma::vec mu, arma::mat sigma){
  unsigned d = x.n_elem - 1;
  arma::vec integral = arma::zeros(d);
  arma::mat N_dirichlet_lebesgue = c_dirichlet_alr_approx(x+1);

  // SIGMA_logistic = 1/(2 * pi * prod(composition(0))^2)

  arma::mat N_logistic = c_logistic_alr_approximation(d);
  arma::mat N_dirichlet_aitchison = c_gaussian_division(N_logistic, N_dirichlet_lebesgue);

  //Rcpp::Rcout << N_logistic << std::endl;
  //Rcpp::Rcout << N_dirichlet_lebesgue << std::endl;
  //Rcpp::Rcout << N_dirichlet_aitchison << std::endl;
  arma::mat N_prior(d,d+1);
  N_prior.head_cols(d) = sigma;
  N_prior.col(d) = mu;

  return(c_gaussian_product(N_prior, N_dirichlet_aitchison));
}

//' @export
// [[Rcpp::export]]
arma::vec c_m1_lrnm_gaussian_approx(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv){
  unsigned d = x.n_elem - 1;
  arma::vec integral = arma::zeros(d);


  return integral;
}

//' @export
// [[Rcpp::export]]
arma::mat c_m2_lrnm_gaussian_approx(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv){
  unsigned d = x.n_elem - 1;
  arma::mat integral = arma::zeros(d,d);

  return integral;
}

