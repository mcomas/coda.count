// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>
#include "coda_base.h"
#include "lrnm_utils.h"
#include "lrnm_gaussian_approx.h"

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix c_rmultinomial_Rcpp(NumericMatrix P, NumericVector vsize) {

  int n = P.nrow();
  int k = P.ncol();

  IntegerMatrix result(n, k);
  for(int i=0; i < n; i++){
    NumericVector prob = P(i, _);
    IntegerVector res(k);
    R::rmultinom(vsize(i), prob.begin(), k, res.begin());
    result(i,_) = res;
  }
  return(result);
}

// [[Rcpp::export]]
arma::mat c_rmultinomial(arma::mat P, arma::vec vsize) {

  int n = P.n_rows;
  int k = P.n_cols;

  arma::mat result(n, k);
  for(int i=0; i < n; i++){
    arma::rowvec prob = P.row(i);
    IntegerVector res(k);
    R::rmultinom(vsize(i), prob.memptr(), k, res.begin());
    result.row(i) = as<arma::rowvec>(res);
  }
  return(result);
}

// [[Rcpp::export]]
arma::mat c_rdirichlet(int n, arma::vec alpha){
  int k = alpha.size();
  arma::mat X(n, k);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < k; j++){
      X(i,j) = R::rgamma(alpha(j), 1);
    }
  }
  arma::mat sumX = arma::sum(X, 1);
  for(int i = 0; i < k; i++){
    X.col(i) = X.col(i) / sumX;
  }
  return(X);
}

// [[Rcpp::export]]
List c_rdirichletmultinomial(arma::vec alpha, arma::vec size){
  int n = size.n_elem;

  arma::mat P = c_rdirichlet(n, alpha);
  return(List::create(P, c_rmultinomial(P, size)));
}


// [[Rcpp::export]]
arma::mat c_rnormal(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


// [[Rcpp::export]]
arma::mat c_rnormalSimplex(int n, arma::vec mu, arma::mat sigma, arma::mat Binv) {
  return inv_B_coordinates(c_rnormal(n, mu, sigma), Binv);
}


// [[Rcpp::export]]
List c_rnormalmultinomial(arma::vec mu, arma::mat sigma, arma::vec size, arma::mat Binv){
  int n = size.n_elem;

  arma::mat P = inv_B_coordinates(c_rnormal(n, mu, sigma), Binv);
  return(List::create(P, c_rmultinomial(P, size)));
}


//' @export
// [[Rcpp::export]]
arma::mat c_rlrnm_posterior(int n, arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv){
  // l_lrnm_join_vec(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv)
  // arma::mat c_posterior_approximation_vec(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){
  int d = x.n_elem - 1;
  arma::mat inv_sigma = arma::inv_sympd(sigma);
  arma::mat N_approx = c_posterior_approximation_vec(x, mu, inv_sigma, Binv);
  arma::mat sigma_proposal = N_approx.head_cols(d);
  arma::vec current = N_approx.col(d);
  double P_current = exp(l_lrnm_join_vec(current, x, mu, inv_sigma, Binv));
  double P_next = 0;
  arma::mat sample = arma::mat(d, n);
  int acceptance = 0;
  for(int i = 0; i < n; i++){
    arma::mat rN = c_rnormal(100, arma::zeros(d), sigma_proposal).t();
    arma::vec unif = arma::randu(100);

    for(int j = 0; j < 100; j++){
      arma::vec next = current + rN.col(j);
      P_next = exp(l_lrnm_join_vec(next, x, mu, inv_sigma, Binv));
      if(unif(j) <= P_next/P_current){
        acceptance++;
        current = arma::vec(next);
        P_current = P_next;
      }
    }
    sample.col(i) = arma::vec(current);
  }
  Rcpp::Rcout << "Acceptance rate:" << ((double)acceptance)/(100*n) << std::endl;
  return(sample.t());
}
