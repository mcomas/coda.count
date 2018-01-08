#include <RcppArmadillo.h>

double dnm(arma::vec x, arma::vec mu, arma::mat sigma, unsigned int order);

arma::vec m1_dnm(arma::vec x, arma::vec mu, arma::mat sigma, unsigned int order);

arma::mat m2_dnm(arma::vec x, arma::vec mu, arma::mat sigma, unsigned int order);

