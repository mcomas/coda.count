// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

//' @export
// [[Rcpp::export]]
arma::vec ldnormal(arma::mat H, arma::vec mu, arma::mat inv_sigma){
  int k = H.n_cols;
  int n = H.n_rows;
  double log_det_val;
  double sign;
  log_det(log_det_val, sign, inv_sigma);
  double norm_const = -0.5 * k * log(2*PI) + 0.5 * log_det_val;
  arma::vec log_norm = arma::vec(n);

  for(int i = 0; i < n; i++){
    arma::mat cal = (H.row(i) - mu.t()) * inv_sigma * (H.row(i) - mu.t()).t();
    log_norm[i] = -0.5 * cal(0);
  }
  arma::vec norm = norm_const + log_norm;

  return(norm);
}

//' @export
// [[Rcpp::export]]
double lpmultinomial_const(arma::vec x){
  int K = x.size();
  double x_total = 0;
  for(int i = 0; i < K; i++) x_total += x[i];
  arma::vec log_sums = arma::vec(1+x_total);
  log_sums[0] = 0;
  for(int l = 1; l <= x_total; l++) log_sums[l] = log_sums[l-1] + log(l);
  double constant = log_sums[x_total];
  for(int i = 0; i < K; i++) constant -= log_sums[x[i]];
  return(constant);
}

double lpmultinomial(arma::vec x, arma::vec p, double lconst){
  //double lconst = lpmultinomial_const(x);
  return( lconst + arma::dot(log(p),x) );
}

//' @export
// [[Rcpp::export]]
arma::mat lpnm_join(arma::vec x, arma::vec mu, arma::mat inv_sigma, arma::mat P, arma::mat H){
  double lconst = lpmultinomial_const(x);
  arma::mat lmult = arma::sum( arma::repmat(arma::mat(x).t(), P.n_rows, 1) % log(P), 1);
  arma::mat lnormal = ldnormal(H, mu, inv_sigma);
  return(lconst + lmult + lnormal);
}

//' @export
// [[Rcpp::export]]
double lpnm_join_no_constant(arma::vec x, arma::vec mu, arma::mat inv_sigma,
                             arma::vec p, arma::vec h){
  double lmult = arma::accu(x % log(p));
  arma::vec y = h-mu;
  double lnormal = -0.5 * ((arma::mat)(y.t() * inv_sigma * y))(0,0);
  return(lmult + lnormal);
}
