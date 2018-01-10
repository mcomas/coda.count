// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
# include "nm.h"
# include "nm_expect.h"
# include "coda.h"
#include "hermite.h"
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

arma::vec ldnormal(arma::mat H, arma::vec mu, arma::mat inv_sigma){
  int k = H.n_cols;
  int n = H.n_rows;

  double norm_const = -0.5 * k * log(2*PI) + 0.5 * log(det(inv_sigma));
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

arma::mat lpnm_join(arma::vec x, arma::vec mu, arma::mat inv_sigma, arma::mat P, arma::mat H){
  double lconst = lpmultinomial_const(x);
  arma::mat lmult = arma::sum( arma::repmat(arma::mat(x).t(), P.n_rows, 1) % log(P), 1);
  arma::mat lnormal = ldnormal(H, mu, inv_sigma);
  return(lconst + lmult + lnormal);
}

arma::vec lpmultinomial_mult(arma::mat P, arma::vec x){
  return(arma::sum( arma::repmat(arma::mat(x).t(), P.n_rows, 1) % log(P), 1) );
}

double lpnm_join_no_constant(arma::vec x, arma::vec mu, arma::mat inv_sigma,
                             arma::vec p, arma::vec h){
  double lmult = arma::accu(x % log(p));
  arma::vec y = h-mu;
  double lnormal = -0.5 * ((arma::mat)(y.t() * inv_sigma * y))(0,0);
  return(lmult + lnormal);
}

// [[Rcpp::export]]
double c_dlrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma,
                       int order = 100, int step_by = 100,
                       double eps = 0.000001, int max_steps = 10){
  double vcurrent = -1;
  double vnext = c_dnm_hermite(x, mu, sigma, order);
  int step = 1;
  while(abs(vcurrent - vnext) > eps &  step < max_steps){
    step++;
    order+=step_by;
    vcurrent = vnext;
    vnext = c_dnm_hermite(x, mu, sigma, order);
  }
  return(vnext);
}

//' @export
// [[Rcpp::export]]
arma::vec c_m1_hermite(arma::vec x, arma::vec mu, arma::mat sigma,
                       int order = 100, int step_by = 100,
                       double eps = 0.000001, int max_steps = 10){
  arma::vec vcurrent(mu.n_elem);
  arma::vec vnext = m1_lrnm_hermite(x, mu, sigma, order);
  int step = 1;
  while(max(abs(vcurrent - vnext)) > eps &  step < max_steps){
    step++;
    order+=step_by;
    vcurrent = vnext;
    vnext = m1_lrnm_hermite(x, mu, sigma, order);
  }
  return(vnext);
}

//' @export
// [[Rcpp::export]]
arma::mat c_m2_hermite(arma::vec x, arma::vec mu, arma::mat sigma,
                       int order = 100, int step_by = 100,
                       double eps = 0.000001, int max_steps = 10){
  arma::mat vcurrent(mu.n_elem, mu.n_elem);
  arma::mat vnext = m2_lrnm_hermite(x, mu, sigma, order);
  int step = 1;
  while(max(max(abs(vcurrent - vnext))) > eps &  step < max_steps){
    step++;
    order+=step_by;
    vcurrent = vnext;
    vnext = m2_lrnm_hermite(x, mu, sigma, order);
  }
  return(vnext);
}

//' @export
// [[Rcpp::export]]
Rcpp::List c_lrnm_fit_hermite(arma::mat X, arma::vec mu0, arma::mat sigma0,
                   int order = 100, int step_by = 100,
                   double eps = 0.000001, int max_steps = 10,
                   int em_max_steps = 10){
  X = X.t();
  int n = X.n_cols;

  arma::mat H(mu0.n_elem, n);
  arma::vec M1(mu0.n_elem);
  arma::mat M2(mu0.n_elem, mu0.n_elem);
  arma::vec mu = mu0;
  arma::mat sigma = sigma0;
  arma::vec mu_prev = mu0 + 1;


  int step = 0;
  while(max(abs(mu_prev - mu)) > eps & step < em_max_steps){
    step++;
    M1.zeros();
    M2.zeros();
    for(int i = 0; i < n; i++){
      double prob = c_dlrnm_hermite(X.col(i), mu, sigma, order, step_by, eps, max_steps);
      H.col(i) = c_m1_hermite(X.col(i), mu, sigma, order, step_by, eps, max_steps)/prob;
      M1 += H.col(i);
      M2 += (c_m2_hermite(X.col(i), mu, sigma, order, step_by, eps, max_steps)/prob);
    }
    mu_prev = mu;
    mu = M1 / n;
    sigma = M2/n - mu * mu.t();
  }

  return Rcpp::List::create(mu, sigma, inv_ilr_coordinates(H.t()));
}

// //' @export
// // [[Rcpp::export]]
// Rcpp::List c_nm_fit_1(arma::mat X, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z,
//                       arma::mat mu_exp, arma::vec var_exp){
//
//   int n = X.n_rows;
//   for(int i = 0; i++; i < n){
//     E = expected_montecarlo_01(X.row(i), mu_ilr, sigma_ilr, Z, mu_exp = H.row(i));
//   }
//   return Rcpp::List::create(alpha, iter);
// }
