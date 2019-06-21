// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "hermite.h"
#include "lrnm_utils.h"
#include "coda_base.h"

//' @export
// [[Rcpp::export]]
double c_d_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv,
                        int order){

  unsigned d = x.n_elem - 1;
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

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
    arma::vec p = exp(Binv * h);
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

//' @export
// [[Rcpp::export]]
arma::vec c_m1_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma,
                            arma::mat Binv, int order){
  unsigned d = x.n_elem - 1;
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;
  arma::vec integral = arma::zeros(d);
  double l_cmult = lpmultinomial_const(x);

  do{
    double w = 0;
    arma::vec h(d);
    for(unsigned int i = 0; i < d; i++){
      h(i) = uni_hermite(index[i],0);
      w += uni_hermite(index[i],1);
    }
    h = mu + rotation * h;
    arma::vec p = exp(Binv * h);
    integral += h * exp(w + lpmultinomial(x, p/accu(p), l_cmult));

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

//' @export
// [[Rcpp::export]]
arma::mat c_m2_lrnm_hermite(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Binv,
                            int order){
  unsigned d = x.n_elem - 1;
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, sigma);
  arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

  unsigned int index[d+1];
  for(unsigned int i = 0; i <= d; i++) index[i] = 0;
  int position = 0, k = 0;
  arma::mat integral = arma::zeros(d,d);
  double l_cmult = lpmultinomial_const(x);

  do{
    double w = 0;
    arma::vec h(d);
    for(unsigned int i = 0; i < d; i++){
      h(i) = uni_hermite(index[i],0);
      w += uni_hermite(index[i],1);
    }
    h = mu + rotation * h;
    arma::vec p = exp(Binv * h);
    integral += h * h.t()* exp(w + lpmultinomial(x, p/accu(p), l_cmult));

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


/*
 * Maybe following functions should be removed
 */

// [[Rcpp::export]]
double c_d_lrnm_hermite_(arma::vec x, arma::vec mu, arma::mat sigma,
                         arma::mat Binv, int order = 100, int step_by = 100,
                         double eps = 0.000001, int max_steps = 10){
  double vcurrent = -1;
  double vnext = c_d_lrnm_hermite(x, mu, sigma, Binv, order);
  int step = 1;
  while(abs(vcurrent - vnext) > eps &  step < max_steps){
    step++;
    order+=step_by;
    vcurrent = vnext;
    vnext = c_d_lrnm_hermite(x, mu, sigma, Binv, order);
  }
  return(vnext);
}

// [[Rcpp::export]]
arma::vec c_m1_lrnm_hermite_(arma::vec x, arma::vec mu, arma::mat sigma,
                             arma::mat Binv, int order = 100, int step_by = 100,
                             double eps = 0.000001, int max_steps = 10){
  arma::vec vcurrent(mu.n_elem);
  arma::vec vnext = c_m1_lrnm_hermite(x, mu, sigma, Binv, order);
  int step = 1;
  while(max(abs(vcurrent - vnext)) > eps &  step < max_steps){
    step++;
    order+=step_by;
    vcurrent = vnext;
    vnext = c_m1_lrnm_hermite(x, mu, sigma, Binv, order);
  }
  return(vnext);
}

// [[Rcpp::export]]
arma::mat c_m2_lrnm_hermite_(arma::vec x, arma::vec mu, arma::mat sigma,
                             arma::mat Binv, int order = 100, int step_by = 100,
                             double eps = 0.000001, int max_steps = 10){
  arma::mat vcurrent(mu.n_elem, mu.n_elem);
  arma::mat vnext = c_m2_lrnm_hermite(x, mu, sigma, Binv, order);
  int step = 1;
  while(max(max(abs(vcurrent - vnext))) > eps &  step < max_steps){
    step++;
    order+=step_by;
    vcurrent = vnext;
    vnext = c_m2_lrnm_hermite(x, mu, sigma, Binv, order);
  }
  return(vnext);
}


// [[Rcpp::export]]
Rcpp::List c_fit_lrnm_hermite_(arma::mat X, arma::vec mu0, arma::mat sigma0,
                               arma::mat Binv, int order = 100, int step_by = 100,
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
      double prob = c_d_lrnm_hermite_(X.col(i), mu, sigma, Binv, order, step_by, eps, max_steps);
      H.col(i) = c_m1_lrnm_hermite_(X.col(i), mu, sigma, Binv, order, step_by, eps, max_steps)/prob;
      M1 += H.col(i);
      M2 += (c_m2_lrnm_hermite_(X.col(i), mu, sigma, Binv, order, step_by, eps, max_steps)/prob);
    }
    mu_prev = mu;
    mu = M1 / n;
    sigma = M2/n - mu * mu.t();
  }

  return Rcpp::List::create(mu, sigma, H.t());
}
