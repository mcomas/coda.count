#include <RcppArmadillo.h>

arma::vec mvf_maximum(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat B,
                      double eps, int max_iter, double prop);

arma::vec mvf_maximum_alr(arma::vec x, arma::vec mu_alr, arma::mat& inv_sigma_alr,
                          arma::vec a, double eps, int max_iter);
