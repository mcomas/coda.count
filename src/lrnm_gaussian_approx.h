#include <RcppArmadillo.h>

arma::mat c_posterior_approximation(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv);
