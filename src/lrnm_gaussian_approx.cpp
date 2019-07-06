// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "gaussian_approx.h"
#include "lrnm_utils.h"
#include "coda_base.h"
#include "dm.h"

arma::mat c_posterior_approximation_vec(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){
  unsigned d = x.n_elem - 1;

  arma::mat N_posterior(d,d+1);
  N_posterior.col(d) = l_lrnm_join_maximum(x, mu, inv_sigma, Binv, 1e-5, 1000);
  //Rcpp::Rcout << N_posterior.col(d);
  arma::mat D2 = l_lrnm_join_d2(N_posterior.col(d), x, mu, inv_sigma, Binv);
  //Rcpp::Rcout << D2;
  N_posterior.head_cols(d) = arma::inv_sympd(-D2);

  return(N_posterior);
}

// [[Rcpp::export]]
arma::cube c_posterior_approximation(arma::mat X, arma::vec mu, arma::mat &sigma, arma::mat &B){
  arma::mat Binv = pinv(B).t();
  arma::mat inv_sigma = inv_sympd(sigma);
  arma::mat Xt = X.t();
  arma::cube approx = arma::cube(mu.n_elem, X.n_cols, X.n_rows);
  for(int i = 0; i < X.n_rows; i++){
    approx.slice(i) = c_posterior_approximation_vec(Xt.col(i), mu, inv_sigma, Binv);
  }
  return(approx);
}

// [[Rcpp::export]]
Rcpp::List c_fit_lrnm_gaussian_approx(arma::mat X, arma::mat B,
                                      double eps = 0.00001, int max_iter = 200, int em_max_steps = 10){

  int n = X.n_rows;
  int d = X.n_cols - 1;
  arma::vec alpha = c_dm_fit_alpha(X);
  arma::mat Binv = pinv(B).t();
  arma::mat P = arma::mat(X);
  P.each_row() += alpha.t();

  arma::mat H = arma::log(P) * B;
  arma::vec mu = mean(H,0).t();
  arma::vec mu_prev;
  arma::mat sigma = cov(H);
  int current_iter = 0;
  do{
    arma::mat inv_sigma = arma::inv_sympd(sigma);
    current_iter++;

    mu_prev = arma::vec(mu);
    arma::vec M1 = arma::zeros(d);
    arma::mat M2 = arma::zeros(d, d);
    for(int i = 0; i < H.n_rows; i++){
      arma::mat N_posterior = c_posterior_approximation_vec(X.row(i).t(), mu, inv_sigma, Binv);
      M1 += N_posterior.col(d);
      M2 += (N_posterior.head_cols(d) + N_posterior.col(d) * N_posterior.col(d).t());
    }
    mu = M1 / n;
    sigma = M2 / n - mu * mu.t();
  } while ( norm(mu-mu_prev, 2) > eps && current_iter < max_iter);

  return Rcpp::List::create(mu, sigma, current_iter);
}
