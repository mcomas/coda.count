// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "lrnm_utils.h"
#include "coda_base.h"
#include "lrnm_gaussian_approx.h"
#include "dm.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List c_fit_lrnm_lm_laplace(arma::mat Y, arma::mat B, arma::mat X,
                                 double eps, int max_iter, arma::mat H0){

  int n = Y.n_rows;
  int k = X.n_cols;
  int d = Y.n_cols - 1;

  arma::mat Binv = pinv(B).t();

  arma::mat H = H0;

  arma::mat beta = arma::inv(X.t() * X) * X.t() * H, beta_prev;
  arma::mat R = H - X * beta;
  arma::mat sigma_lm = R.t() * R / (n-k), inv_sigma_lm;


  int current_iter = 0;
  arma::vec eigval;
  arma::mat eigvec;
  do{
    Rcpp::checkUserInterrupt();

    current_iter++;

    bool cor = eig_sym(eigval, eigvec, sigma_lm);
    if(eigval.min() < 1e-10){
      Rcpp::Rcout << "Covariance matrix is degenerate" << std::endl;
      return Rcpp::List::create(beta, sigma_lm, H, current_iter, eigval, sigma_lm);
    }
    inv_sigma_lm = arma::inv_sympd(sigma_lm);
    beta_prev = arma::mat(beta);
    arma::vec M1 = arma::zeros(d);
    arma::mat M2 = arma::zeros(d, d);
    for(int i = 0; i < Y.n_rows; i++){
      arma::mat mu_i =  X.row(i) * beta;
      arma::mat moments = c_lrnm_posterior_approximation_vec(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
      H.row(i) = moments.col(d).t();
      M1 += moments.col(d);
      M2 += moments.head_cols(d);
    }
    beta = arma::inv(X.t() * X) * X.t() * H;
    R = H - X * beta;
    sigma_lm = R.t() * R / (n-k); // M2 / n;
    //sigma = M2 / n - mu * mu.t();
  } while ( norm(beta-beta_prev, 1) > eps && current_iter < max_iter);

  // Last iteration
  inv_sigma_lm = arma::inv_sympd(sigma_lm);
  for(int i = 0; i < Y.n_rows; i++){
    arma::mat mu_i =  X.row(i) * beta;
    arma::mat moments = c_lrnm_posterior_approximation_vec(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
    H.row(i) = moments.col(d).t();
  }
  beta = arma::inv(X.t() * X) * X.t() * H;

  return Rcpp::List::create(beta, sigma_lm, H, current_iter);
}
