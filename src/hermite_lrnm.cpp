// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "hermite.h"
#include "lrnm_utils.h"
#include "coda_base.h"

// [[Rcpp::export]]
double c_dlrnm_hermite_(arma::vec x, arma::vec mu, arma::mat sigma, int order){
  unsigned d = x.n_elem - 1;
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));
  arma::mat B = ilr_basis(d+1);

  // arma::vec B_mu = B * mu;
  // arma::mat B_rotation = B * rotation;

  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;
  double integral = 0;
  double l_cmult = lpmultinomial_const(x);

  do{
    double w = 0;
    arma::vec h(d);
    for(unsigned int i = 0; i < d; i++){
      h(i) = uni_hermite(index[i],0);
      w += uni_hermite(index[i],1);
    }
    h = mu + rotation * h;
    //arma::vec p = exp(B * (mu + rotation * h));
    arma::vec p = exp(B * h);
    integral += exp(w + lpmultinomial(x, p/accu(p), l_cmult));
    // Calculate next coordinate
    index[position]++;
    while(index[position] == order){
      index[position] = 0;
      position++;
      index[position]++;
    }
    position = 0;
    k++;
  } while (index[d] == 0);

  return integral;
}

// [[Rcpp::export]]
double c_dlrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma,
                       int order = 100, int step_by = 100,
                       double eps = 0.000001, int max_steps = 10){
  double vcurrent = -1;
  double vnext = c_dlrnm_hermite_(x, mu, sigma, order);
  int step = 1;
  while(abs(vcurrent - vnext) > eps &  step < max_steps){
    step++;
    order+=step_by;
    vcurrent = vnext;
    vnext = c_dlrnm_hermite_(x, mu, sigma, order);
  }
  return(vnext);
}
