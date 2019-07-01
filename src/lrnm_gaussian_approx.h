#include <RcppArmadillo.h>

#ifndef LRNM_GAUSSIAN_APPROX_H
#define LRNM_GAUSSIAN_APPROX_H

arma::mat c_posterior_approximation(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv);

#endif
