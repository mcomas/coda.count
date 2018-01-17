// [[Rcpp::depends(RcppArmadillo)]]

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
* - expected_hermite: exact method using hermite poylinomials.
* - expected_montecarlo: montecarlo approximation centered at mu_sampling with variance sigma_sampling.
* - expected_metropolis: montecarlo approximation centered at mu_exp. sigma_exp as sampling variance.
*/

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
Rcpp::List expected_hermite(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, int order){
  double M0 = c_dnm_hermite(x, mu_ilr, sigma_ilr, order);
  arma::vec M1 = m1_lrnm_hermite(x, mu_ilr, sigma_ilr, order)/M0;
  arma::mat M2 = m2_lrnm_hermite(x, mu_ilr, sigma_ilr, order)/M0;
  return(Rcpp::List::create(M0, M1, M2));
}

//' @export
// [[Rcpp::export]]
Rcpp::List expected_montecarlo(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat& Z,
                               arma::vec mu_sampling, arma::mat sigma_sampling, arma::mat& Hz){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::mat B = ilr_basis(K);

  arma::mat inv_sigma_ilr = arma::mat(k,k);
  bool b_inv_sigma_ilr = inv_sympd(inv_sigma_ilr, sigma_ilr);
  if(b_inv_sigma_ilr == false){
    inv_sigma_ilr= inv(diagmat(abs(diagvec(sigma_ilr))));
  }

  arma::mat inv_sigma_sampling = arma::mat(k,k);
  bool b_inv_sigma_sampling = inv_sympd(inv_sigma_sampling, sigma_sampling);
  if(b_inv_sigma_sampling == false){
    sigma_sampling = diagmat(abs(diagvec(sigma_sampling)));
    inv_sigma_sampling = inv(sigma_sampling);
  }

  Hz = Z * arma::chol(sigma_sampling);
  Hz.each_row() += mu_sampling.t();
  arma::vec loglik = -ldnormal(Hz, mu_sampling, inv_sigma_sampling);

  for(int i = 0; i < nsim; i++){
    arma::vec h = Hz.row(i).t();
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
    arma::vec C = Hz.col(i) % lik_st;
    M1(i) += mean(C);
    for(int j = 0;j < k; j++){
      M2(i,j) += mean(C % Hz.col(j));
    }
  }

  return(Rcpp::List::create(M0, M1, M2));

}

//' @export
// [[Rcpp::export]]
Rcpp::List expected_metropolis(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::vec mu_exp,
                               int nsim, int ignored_steps = 100){
  int maxZ = 10000;
  int maxU = 10000;

  int K = x.size();
  int k = K - 1;

  arma::mat Z = arma::randn(k, maxZ);
  arma::vec U = arma::randu(maxU);

  arma::mat inv_sigma = inv_sympd(sigma_ilr);
  arma::mat B = ilr_basis(K);

  // initialisation
  arma::vec h = mu_exp, h_proposal;
  arma::vec p = exp(B * h), p_proposal;

  double f = lpnm_join_no_constant(x, mu_ilr, inv_sigma, p / arma::accu(p), h);

  int irand_z = 0, irand_u = 0, step = 0, nsim_real = 0;
  bool repeat = true;
  double M0 = 0;
  arma::vec M1 = arma::zeros(k);
  arma::mat M2 = arma::mat(k,k);

  for(int i = 0; i < nsim; i++){
    for(int j = 0; j < ignored_steps; j++){
      h_proposal = h + Z.col(irand_z++);
      p_proposal = exp(B * h_proposal);
      double f_proposal = lpnm_join_no_constant(x, mu_ilr, inv_sigma,
                                                p_proposal / arma::accu(p_proposal), h_proposal);

      double cmean = 0.5 * f_proposal + 0.5 * f;
      double alpha = exp(f_proposal-cmean) / exp(f-cmean);
      if(1 < alpha){
        f = f_proposal;
        h = h_proposal;
      }else{
        if(U(irand_u++) < alpha){
          f = f_proposal;
          h = h_proposal;
        }else{
          f = f;
          h = h;
        }
      }
      if(irand_z == maxZ){
        Z = arma::randn(k, maxZ);
        irand_z = 0;
      }
      if(irand_u == maxU){
        U = arma::randu(maxU);
        irand_u = 0;
      }
    }

    M0 += f;
    M1 += h;
    M2 += h * h.t();
  }
  M0 /= nsim;
  M1 /= nsim;
  M2 /= nsim;

  return(Rcpp::List::create(M0, M1, M2));
}

//' @export
// [[Rcpp::export]]
arma::mat metropolis_sample(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::vec x0,
                            int nsim, int ignored_steps = 100){
  int maxZ = 10000;
  int maxU = 10000;

  int K = x.size();
  int k = K - 1;

  arma::mat Z = arma::randn(k, maxZ);
  arma::vec U = arma::randu(maxU);

  arma::mat inv_sigma = inv_sympd(sigma_ilr);
  arma::mat B = ilr_basis(K);

  // initialisation
  arma::vec h = x0, h_proposal;
  arma::vec p = exp(B * h), p_proposal;

  double f = lpnm_join_no_constant(x, mu_ilr, inv_sigma, p / arma::accu(p), h);

  int irand_z = 0, irand_u = 0, step = 0, nsim_real = 0;
  bool repeat = true;

  arma::mat X = arma::mat(k, nsim);

  for(int i = 0; i < nsim; i++){
    for(int j = 0; j < ignored_steps; j++){
      h_proposal = h + Z.col(irand_z++);
      p_proposal = exp(B * h_proposal);
      double f_proposal = lpnm_join_no_constant(x, mu_ilr, inv_sigma,
                                                p_proposal / arma::accu(p_proposal), h_proposal);

      double cmean = 0.5 * f_proposal + 0.5 * f;
      double alpha = exp(f_proposal-cmean) / exp(f-cmean);
      if(1 < alpha){
        f = f_proposal;
        h = h_proposal;
      }else{
        if(U(irand_u++) < alpha){
          f = f_proposal;
          h = h_proposal;
        }else{
          f = f;
          h = h;
        }
      }
      if(irand_z == maxZ){
        Z = arma::randn(k, maxZ);
        irand_z = 0;
      }
      if(irand_u == maxU){
        U = arma::randu(maxU);
        irand_u = 0;
      }
    }

    X.col(i) = h;
  }
  return(X.t());
}
