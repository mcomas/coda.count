#include <RcppArmadillo.h>

arma::vec mvf_maximum(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat B, double eps, int max_iter, double prop);
