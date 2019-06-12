#include <RcppArmadillo.h>

double c_dnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, int order);

arma::vec m1_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, unsigned int order);

arma::mat m2_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, unsigned int order);

