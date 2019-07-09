// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "lrnm_utils.h"
#include "coda_base.h"
#include "lrnm_gaussian_approx.h"
#include "dm.h"

using namespace Rcpp;

// [[Rcpp::export]]
double c_d_lrnm_montecarlo(arma::vec x, arma::vec mu_prior, arma::mat sigma_prior, arma::mat Binv,
                           arma::mat &Z){
  unsigned d = x.n_elem - 1;
  unsigned n = Z.n_rows;

  arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior);
  arma::mat N_posterior = c_posterior_approximation_vec(x, mu_prior, inv_sigma_prior, Binv);
  arma::vec mu = N_posterior.col(d);
  arma::mat sigma = N_posterior.head_cols(d);
  arma::mat inv_sigma = arma::inv_sympd(sigma);


  arma::mat Hz = Z * arma::chol(sigma);

  double l_cmult = l_multinomial_const(x);

  arma::vec h, p;
  double integral = 0;
  for(int i=0; i<n; i++){
    h = Hz.row(i).t() + mu;
    p = arma::exp(Binv * h);

    integral += exp(l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
      l_dnormal_vec(h, mu, inv_sigma) +
      l_multinomial(x, p/accu(p), l_cmult));
  }
  return integral/n;
}


// [[Rcpp::export]]
arma::mat c_moments_lrnm_montecarlo(arma::vec x,
                                    arma::vec mu, arma::mat sigma,
                                    arma::vec mu_prior, arma::mat inv_sigma_prior,
                                    arma::mat Binv, arma::mat &Z){
  unsigned d = x.n_elem - 1;
  unsigned n = Z.n_cols;

  arma::mat inv_sigma = arma::inv_sympd(sigma);


  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;

  arma::mat Hz = arma::chol(sigma).t() * Z;

  double l_cmult = l_multinomial_const(x);

  arma::vec h, p;
  double M0 = 0;
  arma::vec M1 = arma::zeros(d);
  arma::mat M2 = arma::zeros(d,d);

  arma::mat INV_SIGMA = inv_sigma_prior - inv_sigma;
  arma::mat MU = inv_sympd(INV_SIGMA) * (inv_sigma_prior * mu_prior - inv_sigma * mu);

  for(int i=0; i<n; i++){
    h = Hz.col(i) + mu;
    p = arma::exp(Binv * h);

    double dens = exp(l_dnormal_prop_vec(h, MU, INV_SIGMA) +
                      arma::dot(log(p/accu(p)),x));

    // double dens = exp(l_dnormal_prop_vec(h, mu_prior, inv_sigma_prior) -
    //                   l_dnormal_prop_vec(h, mu, inv_sigma) +
    //                   arma::dot(log(p/accu(p)),x));

    M0 += dens;
    M1 += h * dens;
    M2 += (h-mu_prior) * (h-mu_prior).t() * dens;
  }

  arma::mat moments(d, d+1);
  moments.col(d) = M1/M0;
  moments.head_cols(d) = M2/M0;
  return moments;
}

// //' @export
// // [[Rcpp::export]]
// Rcpp::List c_obtain_moments_lrnm_montecarlo(arma::mat Y,
//                                             arma::vec mu, arma::mat sigma,
//                                             arma::mat B, arma::mat &Z){
//   int n = Y.n_rows;
//   int d = Y.n_cols - 1;
//
//   arma::mat Binv = pinv(B).t();
//   arma::mat inv_sigma = arma::inv_sympd(sigma);
//
//   arma::mat M1 = arma::zeros(d, n);
//   arma::cube M2 = arma::zeros(d, d, n);
//
//   for(int i = 0; i < Y.n_rows; i++){
//     arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu, inv_sigma, Binv);
//     arma::mat moments = c_moments_lrnm_montecarlo_precision_lm(Y.row(i).t(),
//                                                                N_posterior.col(d), N_posterior.head_cols(d),
//                                                                mu, inv_sigma,
//                                                                Binv, Z);
//     M1.col(i) = moments.col(d);
//     M2.slice(i) = moments.head_cols(d);
//   }
//   return Rcpp::List::create(M1, M2);
// }


// [[Rcpp::export]]
Rcpp::List c_fit_lrnm_lm_montecarlo(arma::mat Y, arma::mat B, arma::mat X, arma::mat &Z,
                                    double eps, int max_iter){

  int n = Y.n_rows;
  int k = X.n_cols;
  int d = Y.n_cols - 1;

  arma::mat Zt = Z.t();

  arma::vec alpha = c_dm_fit_alpha(Y);
  arma::mat Binv = pinv(B).t();
  arma::mat P = arma::mat(Y);
  P.each_row() += alpha.t();

  arma::mat H = arma::log(P) * B;

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
      arma::mat N_posterior = c_posterior_approximation_vec(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
      arma::mat moments = c_moments_lrnm_montecarlo(Y.row(i).t(),
                                                    N_posterior.col(d), N_posterior.head_cols(d),
                                                    mu_i.t(), inv_sigma_lm,
                                                    Binv, Zt);
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
    arma::mat moments = c_moments_lrnm_montecarlo(Y.row(i).t(),
                                                  N_posterior.col(d), N_posterior.head_cols(d),
                                                  mu_i.t(), inv_sigma_lm,
                                                  Binv, Zt);
    H.row(i) = moments.col(d).t();
  }
  beta = arma::inv(X.t() * X) * X.t() * H;

  return Rcpp::List::create(beta, sigma_lm, H, current_iter);
}