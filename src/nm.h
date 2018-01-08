#include <RcppArmadillo.h>

arma::vec ldnormal(arma::mat H, arma::vec mu, arma::mat inv_sigma);

double lpmultinomial_const(arma::vec x);

double lpmultinomial(arma::vec x, arma::vec p, double lconst);

arma::vec lpmultinomial_mult(arma::mat P, arma::vec x);

double lpnm_join_no_constant(arma::vec x, arma::vec mu, arma::mat inv_sigma, arma::vec p, arma::vec h);
