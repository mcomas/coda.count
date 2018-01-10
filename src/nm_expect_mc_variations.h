#include <RcppArmadillo.h>

arma::vec expected_mc_01_init(arma::vec& x, arma::vec& mu_ilr, arma::mat& sigma_ilr, arma::mat& Z,
                              arma::vec& sampling_mu, arma::mat& Hs);

arma::vec expected_mc_03_init(arma::vec& x, arma::vec& mu_ilr, arma::mat& sigma_ilr,
                              arma::mat& Z, arma::vec& mu_sampling, arma::mat& sigma_sampling, arma::mat& Hs);

arma::vec expected_mc_mean(arma::vec& x, arma::mat& Hs, arma::vec& lik_st);

arma::vec expected_mc_var(arma::vec& x, arma::vec& mu, arma::mat& Hs, arma::vec& lik_st);
