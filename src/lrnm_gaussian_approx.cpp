// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "gaussian_approx.h"
#include "lrnm_utils.h"
#include "coda_base.h"
#include "dm.h"

// //' @export
// // [[Rcpp::export]]
// double c_d_lrnm_gaussian_approx(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv){
//
//   double integral = 0;
//
//   return integral;
// }

//' @export
// [[Rcpp::export]]
arma::mat c_posterior_approximation(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){
  unsigned d = x.n_elem - 1;

  arma::mat N_posterior(d,d+1);
  N_posterior.col(d) = l_lrnm_join_maximum(x, mu, inv_sigma, Binv, 0.000001, 100);
  //Rcpp::Rcout << N_posterior.col(d);
  arma::mat D2 = l_lrnm_join_d2(N_posterior.col(d), x, mu, inv_sigma, Binv);
  //Rcpp::Rcout << D2;
  N_posterior.head_cols(d) = arma::inv_sympd(-D2);

  return(N_posterior);
}

//' @export
// [[Rcpp::export]]
Rcpp::List c_fit_lrnm_gaussian_approx(arma::mat X, arma::mat B,
                                      double eps = 0.00001, int max_iter = 200, int em_max_steps = 10){

  int n = X.n_rows;
  int d = X.n_cols - 1;
  arma::vec alpha = c_dm_fit_alpha(X);
  arma::mat Binv = pinv(B).t();
  arma::mat P = arma::mat(X);
  P.each_row() += alpha.t();

  arma::mat H = arma::log(P) * B;
  arma::vec mu = mean(H,0).t();
  arma::vec mu_prev;
  arma::mat sigma = cov(H);
  int current_iter = 0;
  do{
    arma::mat inv_sigma = arma::inv_sympd(sigma);
    current_iter++;
    //Rcpp::Rcout << "Current: " << current_iter << std::endl;
    mu_prev = arma::vec(mu);
    arma::vec M1 = arma::zeros(d);
    arma::mat M2 = arma::zeros(d, d);
    for(int i = 0; i < H.n_rows; i++){
      arma::mat N_posterior = c_posterior_approximation(X.row(i).t(), mu, inv_sigma, Binv);
      // Rcpp::Rcout << mu << std::endl;
      // Rcpp::Rcout << sigma << std::endl;

      //Rcpp::Rcout << pX << std::endl;
      M1 += N_posterior.col(d);
      M2 += (N_posterior.head_cols(d) + N_posterior.col(d) * N_posterior.col(d).t());
    }
    mu = M1 / n;
    sigma = M2 / n - mu * mu.t();
  } while ( norm(mu-mu_prev, 2) > eps && current_iter < max_iter);

  return Rcpp::List::create(mu, sigma, current_iter);
}

// //' @export
// // [[Rcpp::export]]
// arma::mat c_posterior_alr_approximation(arma::vec x, arma::vec mu, arma::mat sigma){
//   unsigned d = x.n_elem - 1;
//   arma::vec integral = arma::zeros(d);
//   arma::mat N_dirichlet_lebesgue = c_dirichlet_alr_approx(x+1);
//
//   // SIGMA_logistic = 1/(2 * pi * prod(composition(0))^2)
//
//   arma::mat N_logistic = c_logistic_alr_approximation(d);
//   arma::mat N_dirichlet_aitchison = c_gaussian_division(N_logistic, N_dirichlet_lebesgue);
//
//   //Rcpp::Rcout << N_logistic << std::endl;
//   //Rcpp::Rcout << N_dirichlet_lebesgue << std::endl;
//   //Rcpp::Rcout << N_dirichlet_aitchison << std::endl;
//   arma::mat N_prior(d,d+1);
//   N_prior.head_cols(d) = sigma;
//   N_prior.col(d) = mu;
//
//   return(c_gaussian_product(N_prior, N_dirichlet_aitchison));
// }
//
// //' @export
// // [[Rcpp::export]]
// arma::vec c_m1_lrnm_gaussian_approx(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv){
//   unsigned d = x.n_elem - 1;
//   arma::vec integral = arma::zeros(d);
//
//
//   return integral;
// }
//
// //' @export
// // [[Rcpp::export]]
// arma::mat c_m2_lrnm_gaussian_approx(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv){
//   unsigned d = x.n_elem - 1;
//   arma::mat integral = arma::zeros(d,d);
//
//   return integral;
// }

