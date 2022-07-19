// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppEnsmallen.h>
#include "lrnm_utils.h"
#include "lrnm_gaussian_approx.h"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::mat c_moments_lrnm_montecarlo_sigma_inverse(arma::vec x,
                                                  arma::vec mu, arma::mat inv_sigma,
                                                  arma::vec mu_prior, arma::mat inv_sigma_prior,
                                                  arma::mat Binv, arma::mat &Z,
                                                  arma::vec mu_centering){
  unsigned d = x.n_elem - 1;
  unsigned n = Z.n_cols;

  arma::mat Hz = chol(arma::inv_sympd(inv_sigma)).t() * Z;

  arma::vec h, p;
  double M0 = 0;
  arma::vec M1 = arma::zeros(d);
  arma::mat M2 = arma::zeros(d,d);

  arma::mat p_mu = arma::exp(Binv * mu);
  double constant = -arma::dot(log(p_mu/accu(p_mu)),x);
  for(int i=0; i<n; i++){
    h = Hz.col(i) + mu;
    p = arma::exp(Binv * h);
    double dens = exp(l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
                      l_dnormal_vec(h, mu, inv_sigma) +
                      arma::dot(log(p/accu(p)),x) + constant);

    M0 += dens;
    M1 += h * dens;
    M2 += (h-mu_centering) * (h-mu_centering).t() * dens;
  }

  arma::mat moments(d, d+1);
  moments.col(d) = M1/M0;
  moments.head_cols(d) = M2/M0;
  return moments;
}


//' @export
// [[Rcpp::export]]
Rcpp::List c_fit_lrnm_lm_montecarlo(arma::mat Y, arma::mat B, arma::mat X, arma::mat &Z,
                                    double eps, int max_iter, arma::mat H0){

  int n = Y.n_rows;
  int k = X.n_cols;
  int d = Y.n_cols - 1;

  arma::mat Zt = Z.t();


  arma::mat Binv = pinv(B).t();

  arma::mat H = H0;
  arma::mat invXtXX = arma::inv_sympd(X.t() * X) * X.t();
  arma::mat beta = invXtXX * H, beta_prev;
  arma::mat R = H - X * beta;
  arma::mat sigma_lm = R.t() * R / (n-k), inv_sigma_lm;


  int current_iter = 0;
  arma::vec eigval;
  arma::mat eigvec;
  do{
    Rcpp::checkUserInterrupt();

    current_iter++;

    bool cor = eig_sym(eigval, eigvec, sigma_lm);

    if(eigval.min() < 1e-10 | eigval.max() < 1e-5){
      Rcpp::Rcout << "Covariance matrix is degenerate" << std::endl;
      return Rcpp::List::create(beta, sigma_lm, H, current_iter, eigval, sigma_lm);
    }

    inv_sigma_lm = arma::inv_sympd(sigma_lm);

    beta_prev = arma::mat(beta);
    arma::vec M1 = arma::zeros(d);
    arma::mat M2 = arma::zeros(d, d);

    for(int i = 0; i < Y.n_rows; i++){
      if(i % 50 == 0){
        Rcpp::checkUserInterrupt();
      }

      arma::mat mu_i =  X.row(i) * beta;
      arma::mat N_posterior_si = c_posterior_approximation_vec_sigma_inverse(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);

      arma::mat moments = c_moments_lrnm_montecarlo_sigma_inverse(Y.row(i).t(),
                                                                  N_posterior_si.col(d), N_posterior_si.head_cols(d),
                                                                  mu_i.t(), inv_sigma_lm,
                                                                  Binv, Zt, mu_i.t());
      H.row(i) = moments.col(d).t();
      M1 += moments.col(d);
      M2 += moments.head_cols(d);
    }
    beta = invXtXX * H;
    //R = H - X * beta;
    sigma_lm = M2 / n; //R.t() * R / (n-k);
    // Rcpp::Rcout << sigma_lm;
    Rcpp::Rcout << (beta-beta_prev) / beta_prev;
    Rcpp::Rcout << beta << norm(beta-beta_prev, 2) << " " << norm(beta-beta_prev, 1) << std::endl;
    //sigma = M2 / n - mu * mu.t();
  } while ( norm(beta-beta_prev, 1) > eps && current_iter < max_iter);
  // Rcpp::Rcout << inv_sympd(sigma_lm) << std::endl;
  // Last iteration
  inv_sigma_lm = arma::inv_sympd(sigma_lm);
  for(int i = 0; i < Y.n_rows; i++){
    arma::mat mu_i =  X.row(i) * beta;
    arma::mat N_posterior_si = c_posterior_approximation_vec_sigma_inverse(Y.row(i).t(), mu_i.t(), inv_sigma_lm, Binv);
    arma::mat moments = c_moments_lrnm_montecarlo_sigma_inverse(Y.row(i).t(),
                                                                N_posterior_si.col(d), N_posterior_si.head_cols(d),
                                                                mu_i.t(), inv_sigma_lm,
                                                                Binv, Zt, mu_i.t());
    H.row(i) = moments.col(d).t();
  }
  beta = invXtXX * H;

  return Rcpp::List::create(beta, sigma_lm, H, current_iter);
}
