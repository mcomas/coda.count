#include <RcppArmadillo.h>

#ifndef LRNM_GAUSSIAN_APPROX_H
#define LRNM_GAUSSIAN_APPROX_H

//arma::mat c_posterior_approximation(arma::mat x, arma::vec mu, arma::mat &sigma, arma::mat &B, double eps = 1e-05, int niter = 1000);
arma::mat c_lrnm_posterior_approximation_vec(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv, double eps = 1e-05, int niter = 1000);
arma::mat c_lrnm_posterior_approximation_vec_sigma_inverse(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv, double eps, int niter);
arma::mat c_lrnm_cond_posterior_approximation_vec(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv, double eps, int niter);

#endif
