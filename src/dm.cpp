// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
double c_ddm(arma::vec x, arma::vec alpha){
  double left = (R::lgammafn(1 + arma::sum(x)) + R::lgammafn( arma::sum(alpha) )) - R::lgammafn( arma::sum(x+alpha) );
  double right = 0;
  for(int i=0;i<x.size();i++){
    right += R::lgammafn(x[i] + alpha[i]) - ( R::lgammafn(1 + x[i]) + R::lgammafn(alpha[i]) );
  }
  return( exp(left + right ) );
}

// [[Rcpp::export]]
Rcpp::List c_dm_fit(arma::mat& X, double eps, int maxiter){

  int K = X.n_rows;
  int n = X.n_cols;

  arma::rowvec N = arma::sum(X, 0);

  arma::mat P(K,n);
  for(int i=0;i<n;i++){
    P.col(i) = X.col(i) / arma::accu(X.col(i));
  }
  arma::vec alpha_prev = mean(P,1);
  arma::vec m2 = mean(P % P,1);
  double s = arma::median( (alpha_prev - m2) / (m2 - arma::square(alpha_prev)) );
  alpha_prev *= s;

  arma::vec alpha(K);
  int iter = 0;
  double err = eps + 1;
  while(err > eps & iter < maxiter){
    double alpha_total = arma::accu(alpha_prev);
    double psi_alpha_total = R::digamma(alpha_total);

    for(int j=0; j<K; j++){
      double numerator = 0, denominator = 0;
      for(int i=0; i<n; i++){
        numerator += (R::digamma(X(j,i) + alpha_prev(j)) - R::digamma(alpha_prev(j)));
        denominator += (R::digamma(N(i) + alpha_total) - psi_alpha_total);
      }
      alpha(j) = alpha_prev(j) * numerator / denominator;
    }
    err = arma::norm(alpha - alpha_prev);
    alpha_prev = alpha;

    iter++;
  }
  List ret;
  ret["alpha"] = alpha;
  ret["iter"] = iter;
  return ret;
}

