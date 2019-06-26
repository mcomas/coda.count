#include <RcppArmadillo.h>

#ifndef LRNM_UTILS_H
#define LRNM_UTILS_H

double l_multinomial_const(arma::vec x);
double l_multinomial(arma::vec x, arma::vec p, double lconst);
arma::vec l_lrnm_join_d1(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv);
arma::vec l_lrnm_join_d2(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv);
arma::vec l_lrnm_join_maximum(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv,
                              double eps = 0.000001, int max_iter = 100);

#endif
