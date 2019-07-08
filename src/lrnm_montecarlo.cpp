// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "lrnm_utils.h"
#include "coda_base.h"
#include "lrnm_gaussian_approx.h"
#include "dm.h"

using namespace Rcpp;

// [[Rcpp::export]]
double c_d_lrnm_montecarlo(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv,
                           int order){
  double integral = 0;

  return integral;
}


// [[Rcpp::export]]
arma::mat c_moments_lrnm_montecarlo_precision_lm(arma::vec x,
                                                 arma::vec mu, arma::mat sigma,
                                                 arma::vec mu_prior, arma::mat sigma_prior,
                                                 arma::mat Binv, int order){
  arma::mat moments(d, d+1);
  return moments;
}

//' @export
// [[Rcpp::export]]
Rcpp::List c_obtain_moments_lrnm_montecarlo(arma::mat Y,
                                            arma::vec mu, arma::mat sigma,
                                            arma::mat B, int order){
  int n = Y.n_rows;
  int d = Y.n_cols - 1;

  arma::mat Binv = pinv(B).t();
  arma::mat inv_sigma = arma::inv_sympd(sigma);

  arma::mat M1 = arma::zeros(d, n);
  arma::cube M2 = arma::zeros(d, d, n);

  for(int i = 0; i < Y.n_rows; i++){
    arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu, inv_sigma, Binv);
    arma::mat moments = c_moments_lrnm_montecarlo_precision_lm(Y.row(i).t(),
                                                               N_posterior.col(d), N_posterior.head_cols(d),
                                                               mu, sigma,
                                                               Binv, order);
    M1.col(i) = moments.col(d);
    M2.slice(i) = moments.head_cols(d);
  }
  return Rcpp::List::create(M1, M2);
}

//' @export
// [[Rcpp::export]]
Rcpp::List c_fit_lm_lrnm_montecarlo_centered(arma::mat Y, arma::mat B, arma::mat X, int order,
                                             double eps, int max_iter){

  int n = Y.n_rows;
  int k = X.n_cols;
  int d = Y.n_cols - 1;
  arma::vec alpha = c_dm_fit_alpha(Y);
  arma::mat Binv = pinv(B).t();
  arma::mat P = arma::mat(Y);
  P.each_row() += alpha.t();

  arma::mat H = arma::log(P) * B;

  arma::mat beta = arma::inv(X.t() * X) * X.t() * H, beta_prev;
  arma::mat R = H - X * beta;
  arma::mat sigma_lm = R.t() * R / (n-k), inv_sigma_lm;


  int current_iter = 0;
  do{
    current_iter++;

    inv_sigma_lm = arma::inv_sympd(sigma_lm);
    beta_prev = arma::mat(beta);
    arma::vec M1 = arma::zeros(d);
    arma::mat M2 = arma::zeros(d, d);
    for(int i = 0; i < Y.n_rows; i++){
      arma::mat mu_i =  X.row(i) * beta;
      arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
      arma::mat moments = c_moments_lrnm_montecarlo_precision_lm(Y.row(i).t(),
                                                                 N_posterior.col(d), N_posterior.head_cols(d),
                                                                 mu_i.t(), sigma_lm,
                                                                 Binv, order);
      H.row(i) = moments.col(d).t();
      M1 += moments.col(d);
      M2 += moments.head_cols(d);
    }
    beta = arma::inv(X.t() * X) * X.t() * H;
    //R = H - X * beta;
    sigma_lm = M2 / n; //R.t() * R / (n-k);
    //sigma = M2 / n - mu * mu.t();
  } while ( norm(beta-beta_prev, 2) > eps && current_iter < max_iter);

  // Last iteration
  inv_sigma_lm = arma::inv_sympd(sigma_lm);
  for(int i = 0; i < Y.n_rows; i++){
    arma::mat mu_i =  X.row(i) * beta;
    arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
    arma::mat moments = c_moments_lrnm_montecarlo_precision_lm(Y.row(i).t(),
                                                               N_posterior.col(d), N_posterior.head_cols(d),
                                                               mu_i.t(), sigma_lm,
                                                               Binv, order);

    H.row(i) = moments.col(d).t();
  }
  beta = arma::inv(X.t() * X) * X.t() * H;

  return Rcpp::List::create(beta, sigma_lm, H, current_iter);
}
