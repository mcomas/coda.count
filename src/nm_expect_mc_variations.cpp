// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "coda.h"
#include "nm.h"
#include "dm.h"
#include "hermite.h"
#include <random>
#include <vector>


// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/*
*
* - expected_montecarlo_01: montecarlo approximation centered at mu_exp. sigma_ilr used as sampling variance.
* - expected_montecarlo_02: montecarlo approximation centered at mu_exp. sigma_ilr with total variance sigma_exp as sampling variance.
* - expected_montecarlo_03: metropolis sampling
*/

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
Rcpp::List expected_montecarlo_01(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z,
                      arma::vec mu_exp){ //, double var_exp){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;

  arma::mat inv_sigma = inv_sympd(sigma_ilr);

  arma::vec sampling_mu =  mu_exp;
  arma::mat SAMPLING_MU = arma::repmat(sampling_mu.t(), nsim, 1);

  arma::mat sampling_sigma = sigma_ilr;
  arma::mat sampling_inv_sigma = inv_sigma;
  arma::mat sampling_sigma_chol = arma::chol(sampling_sigma);

  //arma::mat mu12 = 0.5 * (sampling_mu.t() * inv_sigma * sampling_mu - mu_ilr.t() * inv_sigma * mu_ilr);
  arma::mat D = inv_sigma * (mu_ilr-sampling_mu);

  //double mult_const = lpmultinomial_const(x);

  double M0 = 0;
  arma::vec M1 = arma::vec(k);
  arma::mat M2 = arma::mat(k,k);
  M1.zeros();M2.zeros();
  arma::mat Hs = SAMPLING_MU + Z * sampling_sigma_chol;
  arma::mat Ps = inv_ilr_coordinates(Hs);
  arma::vec loglik = lpmultinomial_mult(Ps, x) + Hs * D; //mu12(0,0) +
  double cmax = max(loglik);

  arma::vec lik = exp(loglik - cmax);
  arma::vec lik_st = lik / mean(lik);

  M0 += mean(lik);
  for(int i = 0;i < k; i++){
    arma::vec C = Hs.col(i) % lik_st;
    M1(i) += mean(C);
    for(int j = 0;j < k; j++){
      M2(i,j) += mean(C % Hs.col(j));
    }
  }
  return(Rcpp::List::create(M0, M1, M2));
}

arma::vec expected_montecarlo_01_init(arma::vec x, arma::vec& mu_ilr, arma::mat& sigma_ilr, arma::mat& Z,
                                      arma::vec mu_exp, arma::mat& Hs){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;

  arma::mat inv_sigma = inv_sympd(sigma_ilr);

  arma::vec sampling_mu =  mu_exp;
  arma::mat SAMPLING_MU = arma::repmat(sampling_mu.t(), nsim, 1);

  arma::mat sampling_sigma = sigma_ilr;
  arma::mat sampling_inv_sigma = inv_sigma;
  arma::mat sampling_sigma_chol = arma::chol(sampling_sigma);

  arma::mat D = inv_sigma * (mu_ilr-sampling_mu);

  double M0 = 0;
  arma::vec M1 = arma::vec(k);
  M1.zeros();
  Hs = SAMPLING_MU + Z * sampling_sigma_chol;
  arma::mat Ps = inv_ilr_coordinates(Hs);
  arma::vec loglik = lpmultinomial_mult(Ps, x) + Hs * D; //mu12(0,0) +
  double cmax = max(loglik);

  arma::vec lik = exp(loglik - cmax);
  arma::vec lik_st = lik / mean(lik);

  return(lik_st);
}

//' @export
// [[Rcpp::export]]
arma::vec expected_mc_01_init(arma::vec& x, arma::vec& mu_ilr, arma::mat& sigma_ilr,
                              arma::mat& Z, arma::vec& sampling_mu, arma::mat& Hs){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;

  arma::mat sampling_sigma = sigma_ilr;

  arma::mat inv_sigma = arma::mat(K-1,K-1);

  bool b_inv_sigma = inv_sympd(inv_sigma, sigma_ilr);
  if(b_inv_sigma == false){
    inv_sigma = inv(diagmat(abs(diagvec(sigma_ilr))));
  }


  arma::mat sampling_sigma_chol = arma::mat(K-1,K-1);
  bool b_sampling_sigma_chol =  arma::chol(sampling_sigma_chol, sampling_sigma);
  if(b_sampling_sigma_chol == false){
    arma::vec diagonal = abs(diagvec(sampling_sigma));
    diagonal.transform( [](double val) { return (val + 1e-10); } );
    sampling_sigma_chol = diagmat(sqrt(diagonal));
  }

  arma::mat D = inv_sigma * (mu_ilr-sampling_mu);

  Hs = Z * sampling_sigma_chol;
  Hs.each_row() += sampling_mu.t();
  arma::mat Ps = inv_ilr_coordinates(Hs);
  arma::vec loglik = lpmultinomial_mult(Ps, x) + Hs * D; //mu12(0,0) +
  double cmax = max(loglik);

  arma::vec lik = exp(loglik - cmax);
  arma::vec lik_st = lik / mean(lik);

  return(lik_st);
}

//' @export
// [[Rcpp::export]]
arma::vec expected_mc_03_init(arma::vec& x, arma::vec& mu_ilr, arma::mat& inv_sigma_ilr,
                              arma::mat& Z, arma::vec& mu_sampling, arma::mat& sigma_sampling, arma::mat& Hs){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::mat B = ilr_basis(K);

  //Rcpp::Rcout << "in ";
  arma::mat inv_sigma_sampling = arma::mat(k,k);
  bool b_inv_sigma_sampling = inv_sympd(inv_sigma_sampling, sigma_sampling);
  if(b_inv_sigma_sampling == false){
    sigma_sampling = diagmat(abs(diagvec(sigma_sampling)));
    inv_sigma_sampling = inv(sigma_sampling);
  }

  //Rcpp::Rcout << "out " << std::endl;
  arma::mat chol_decomp = arma::mat(k,k);
  bool b_chol_decomp = arma::chol(chol_decomp, sigma_sampling);
  if(b_chol_decomp == false){
    Rcpp::Rcout << sigma_sampling << std::endl;
    arma::vec diagonal = abs(diagvec(sigma_sampling));
    diagonal.transform( [](double val) { return (val + 1e-10); } );
    chol_decomp = diagmat(sqrt(diagonal));
  }

  Hs = Z * chol_decomp;
  Hs.each_row() += mu_sampling.t();
  arma::vec loglik = -ldnormal(Hs, mu_sampling, inv_sigma_sampling);

  for(int i = 0; i < nsim; i++){
    arma::vec h = Hs.row(i).t();
    arma::vec p = exp(B * h);
    loglik(i) += lpnm_join_no_constant(x, mu_ilr, inv_sigma_ilr, p / arma::accu(p), h);
  }
  arma::vec lik = arma::exp(loglik - max(loglik));
  arma::vec lik_st = lik / mean(lik);


  return(lik_st);

}

//' @export
// [[Rcpp::export]]
arma::vec expected_mc_mean(arma::vec& x, arma::mat& Hs, arma::vec& lik_st){
  int k = x.size() - 1;
  arma::vec M1 = arma::vec(k);
  M1.zeros();
  for(int i = 0;i < k; i++){
    arma::vec C = Hs.col(i) % lik_st;
    M1(i) += mean(C);
  }
  return(M1);
}

//' @export
// [[Rcpp::export]]
arma::mat expected_mc_var(arma::vec& x, arma::vec& mu, arma::mat& Hs, arma::vec& lik_st){
  int k = x.size() - 1;
  arma::mat M2 = arma::mat(k, k);
  M2.zeros();
  for(int i = 0; i < k; i++){
    for(int j = 0; j < k; j++){
      M2(i,j) += mean( (Hs.col(i) - mu(i)) % (Hs.col(j) - mu(j)) % lik_st);
    }
  }
  return(M2);
}

//' @export
// [[Rcpp::export]]
arma::vec expected_montecarlo_01_mean(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z,
                                      arma::vec mu_exp){
  int k = x.size() - 1;
  arma::mat Hs = arma::mat(Z);
  arma::vec lik_st = expected_montecarlo_01_init(x, mu_ilr, sigma_ilr, Z, mu_exp, Hs);
  arma::vec M1 = arma::vec(k);
  M1.zeros();
  for(int i = 0;i < k; i++){
    arma::vec C = Hs.col(i) % lik_st;
    M1(i) += mean(C);
  }
  return(M1);
}

//' @export
// [[Rcpp::export]]
Rcpp::List expected_montecarlo_01_centered(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z,
                               arma::vec mu_exp){ //, double var_exp){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;

  arma::mat inv_sigma = inv_sympd(sigma_ilr);

  arma::vec sampling_mu =  mu_exp;
  arma::mat SAMPLING_MU = arma::repmat(sampling_mu.t(), nsim, 1);

  arma::mat sampling_sigma = sigma_ilr;
  arma::mat sampling_inv_sigma = inv_sigma;
  arma::mat sampling_sigma_chol = arma::chol(sampling_sigma);

  //arma::mat mu12 = 0.5 * (sampling_mu.t() * inv_sigma * sampling_mu - mu_ilr.t() * inv_sigma * mu_ilr);
  arma::mat D = inv_sigma * (mu_ilr-sampling_mu);

  //double mult_const = lpmultinomial_const(x);

  double M0 = 0;
  arma::vec M1 = arma::vec(k);
  arma::mat M2 = arma::mat(k,k);
  M1.zeros();M2.zeros();
  arma::mat Hs = SAMPLING_MU + Z * sampling_sigma_chol;
  arma::mat Ps = inv_ilr_coordinates(Hs);
  arma::vec loglik = lpmultinomial_mult(Ps, x) + Hs * D; //mu12(0,0) +
  double cmax = max(loglik);

  arma::vec lik = exp(loglik - cmax);
  arma::vec lik_st = lik / mean(lik);

  M0 += mean(lik);
  for(int i = 0;i < k; i++){
    M1(i) += mean( Hs.col(i) % lik_st );
  }
  for(int i = 0;i < k; i++){
    for(int j = 0;j < k; j++){
      M2(i,j) += mean( (Hs.col(i) - M1(i)) % (Hs.col(j) - M1(j)) % lik_st );
    }
  }
  return(Rcpp::List::create(M0, M1, M2));
}

//' @export
// [[Rcpp::export]]
arma::mat expected_montecarlo_02(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z,
                     arma::vec mu_exp, double var_exp){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;

  arma::mat inv_sigma = inv_sympd(sigma_ilr);

  arma::vec sampling_mu =  mu_exp;
  arma::mat SAMPLING_MU = arma::repmat(sampling_mu.t(), nsim, 1);

  arma::mat sampling_sigma = var_exp * sigma_ilr;
  arma::mat sampling_inv_sigma = inv_sigma / var_exp;
  arma::mat sampling_sigma_chol = arma::chol(sampling_sigma);

  //arma::mat mu12 = 0.5 * (sampling_mu.t() * inv_sigma * sampling_mu - mu_ilr.t() * inv_sigma * mu_ilr);
  arma::mat D = inv_sigma * (mu_ilr-sampling_mu);

  //double mult_const = lpmultinomial_const(x);


  arma::mat M = arma::mat(1+k,k);
  M.zeros();
  arma::mat Hs = SAMPLING_MU + Z * sampling_sigma_chol;
  arma::mat Ps = inv_ilr_coordinates(Hs);
  arma::vec loglik = lpmultinomial_mult(Ps, x) +
    0.5 * sum(((1/var_exp-1) * Hs * inv_sigma) % Hs, 1) +
    Hs * inv_sigma * (mu_ilr - 1/var_exp * sampling_mu);
  double cmax = max(loglik);

  arma::vec lik = exp(loglik - cmax);
  arma::vec lik_st = lik / mean(lik);

  for(int i = 0;i < k; i++){
    arma::vec C = Hs.col(i) % lik_st;
    M(0,i) += mean(C);
    for(int j = 0;j < k; j++){
      M(1+i,j) += mean(C % Hs.col(j));
    }
  }
  return(M);
}

//' @export
// [[Rcpp::export]]
Rcpp::List expected_montecarlo_03(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z,
                      arma::vec mu_sampling, arma::mat sigma_sampling){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::mat B = ilr_basis(K);

  arma::mat inv_sigma_ilr = inv_sympd(sigma_ilr);
  arma::mat inv_sigma_sampling = inv_sympd(sigma_sampling);


  Z = Z * arma::chol(sigma_sampling);
  Z.each_row() += mu_sampling.t();
  arma::vec loglik = -ldnormal(Z, mu_sampling, inv_sigma_sampling);

  for(int i = 0; i < nsim; i++){
    arma::vec h = Z.row(i).t();
    arma::vec p = exp(B * h);
    loglik(i) += lpnm_join_no_constant(x, mu_ilr, inv_sigma_ilr, p / arma::accu(p), h);
  }
  arma::vec lik = arma::exp(loglik - max(loglik));
  arma::vec lik_st = lik / mean(lik);

  double M0 = 0;
  arma::vec M1 = arma::vec(k);
  arma::mat M2 = arma::mat(k,k);
  M1.zeros(); M2.zeros();

  M0 += mean(lik);
  for(int i = 0;i < k; i++){
    arma::vec C = Z.col(i) % lik_st;
    M1(i) += mean(C);
    for(int j = 0;j < k; j++){
      M2(i,j) += mean(C % Z.col(j));
    }
  }

  return(Rcpp::List::create(M0, M1, M2));

}

//' @export
// [[Rcpp::export]]
Rcpp::List expected_montecarlo_04(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z,
                                  arma::vec m1, arma::mat m2){
  arma::vec mu_sampling = m1;
  arma::mat sigma_sampling =  m2 - m1 * m1.t();

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::mat B = ilr_basis(K);

  arma::mat inv_sigma_ilr = inv_sympd(sigma_ilr);
  arma::mat inv_sigma_sampling = inv_sympd(sigma_sampling);


  Z = Z * arma::chol(sigma_sampling);
  Z.each_row() += mu_sampling.t();
  arma::vec loglik = -ldnormal(Z, mu_sampling, inv_sigma_sampling);

  for(int i = 0; i < nsim; i++){
    arma::vec h = Z.row(i).t();
    arma::vec p = exp(B * h);
    loglik(i) += lpnm_join_no_constant(x, mu_ilr, inv_sigma_ilr, p / arma::accu(p), h);
  }
  arma::vec lik = arma::exp(loglik - max(loglik));
  arma::vec lik_st = lik / mean(lik);

  double M0 = 0;
  arma::vec M1 = arma::vec(k);
  arma::mat M2 = arma::mat(k,k);
  M1.zeros(); M2.zeros();

  M0 += mean(lik);
  for(int i = 0;i < k; i++){
    arma::vec C = Z.col(i) % lik_st;
    M1(i) += mean(C);
    for(int j = 0;j < k; j++){
      M2(i,j) += mean(C % Z.col(j));
    }
  }

  return(Rcpp::List::create(M0, M1, M2));

}

