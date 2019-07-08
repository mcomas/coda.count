// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "hermite.h"
#include "lrnm_utils.h"
#include "coda_base.h"
#include "lrnm_gaussian_approx.h"
#include "dm.h"

using namespace Rcpp;

// [[Rcpp::export]]
double c_d_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv,
                        int order){

  unsigned d = x.n_elem - 1;
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;
  double integral = 0;
  double l_cmult = l_multinomial_const(x);

  do{
    double w = 0;
    arma::vec h(d);
    for(unsigned int i = 0; i < d; i++){
      h(i) = uni_hermite(index[i],0);
      w += uni_hermite(index[i],1);
    }
    h = mu + rotation * h;
    arma::vec p = exp(Binv * h);
    integral += exp(w + l_multinomial(x, p/accu(p), l_cmult));
    // Calculate next coordinate
    index[position]++;
    while(index[position] == order){
      index[position] = 0;
      position++;
      index[position]++;
    }
    position = 0;
    k++;
  } while (index[d] == 0);

  return integral;
}

// [[Rcpp::export]]
arma::mat c_moments_lrnm_hermite_precision_lm(arma::vec x,
                                              arma::vec mu, arma::mat sigma,
                                              arma::vec mu_prior, arma::mat sigma_prior,
                                              arma::mat Binv, int order){
  unsigned d = x.n_elem - 1;
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  arma::mat inv_sigma = arma::inv_sympd(sigma);
  arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior);

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;
  double M0 = 0;
  arma::vec M1 = arma::zeros(d);
  arma::mat M2 = arma::zeros(d,d);
  double l_cmult = l_multinomial_const(x);

  do{
    double w = 0;
    arma::vec h(d);
    for(unsigned int i = 0; i < d; i++){
      h(i) = uni_hermite(index[i],0);
      w += uni_hermite(index[i],1);
    }
    h = mu + rotation * h;
    arma::vec p = exp(Binv * h);
    double dens = exp(w + ldnormal_vec(h, mu_prior, inv_sigma_prior) -
                      ldnormal_vec(h, mu, inv_sigma) +
                      l_multinomial(x, p/accu(p), l_cmult));
    M0 += dens;
    M1 += h * dens;
    M2 += (h-mu_prior) * (h-mu_prior).t() * dens;

    // Calculate next coordinate
    index[position]++;
    while(index[position] == order){
      index[position] = 0;
      position++;
      index[position]++;
    }
    position = 0;
    k++;
  } while (index[d] == 0);

  arma::mat moments(d, d+1);
  moments.col(d) = M1/M0;
  moments.head_cols(d) = M2/M0;
  return moments;
}

//' @export
// [[Rcpp::export]]
Rcpp::List c_obtain_moments_lrnm_hermite(arma::mat Y,
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
    arma::mat moments = c_moments_lrnm_hermite_precision_lm(Y.row(i).t(),
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
Rcpp::List c_fit_lm_lrnm_hermite_centered(arma::mat Y, arma::mat B, arma::mat X, int order,
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
      arma::mat moments = c_moments_lrnm_hermite_precision_lm(Y.row(i).t(),
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
    arma::mat moments = c_moments_lrnm_hermite_precision_lm(Y.row(i).t(),
                                                            N_posterior.col(d), N_posterior.head_cols(d),
                                                            mu_i.t(), sigma_lm,
                                                            Binv, order);

    H.row(i) = moments.col(d).t();
  }
  beta = arma::inv(X.t() * X) * X.t() * H;

  return Rcpp::List::create(beta, sigma_lm, H, current_iter);
}
