#include <RcppArmadillo.h>

#ifndef LRNM_UTILS_H
#define LRNM_UTILS_H

double lpmultinomial_const(arma::vec x);
double lpmultinomial(arma::vec x, arma::vec p, double lconst);

#endif
