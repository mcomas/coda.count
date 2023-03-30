// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppEnsmallen.h>
#include "dm.h"
#include "lrnm_utils.h"
#include "lrnm_gaussian_approx.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec c_d_lrnm_montecarlo(arma::mat& X,
                              arma::vec& mu, arma::mat& sigma,
                              arma::mat& Binv, arma::mat &Z){
  unsigned d = X.n_rows - 1;
  unsigned n = X.n_cols;
  unsigned nsim = Z.n_cols;
  arma::mat sigma_inv = arma::inv_sympd(sigma);

  arma::vec h, p;
  arma::vec dens = arma::zeros(n);
  for(int i=0;i<n;i++){
    double l_mult_constant;
    l_mult_constant = l_multinomial_const(X.col(i));
    arma::mat N_posterior_si(d,d+1);
    N_posterior_si = c_lrnm_posterior_approximation_vec_sigma_inverse(X.col(i), mu, sigma_inv, Binv, 1e-05, 1000);

    arma::vec mu_laplace = N_posterior_si.col(d);
    arma::mat sigma_inv_laplace = N_posterior_si.head_cols(d);


    arma::mat Hz = chol(arma::inv_sympd(sigma_inv_laplace)).t() * Z;
    for(int k=0;k<nsim;k++){
      h = Hz.col(k) + mu_laplace;
      p = arma::exp(Binv * h);
      dens(i) += exp(l_dnormal_vec(h, mu, sigma_inv) -
        l_dnormal_vec(h, mu_laplace, sigma_inv_laplace) +
        arma::dot(log(p/accu(p)),X.col(i)) + l_mult_constant);
    }
  }
  return(dens/nsim);

}

//' @export
// [[Rcpp::export]]
List c_lrnm_posterior_moments_montecarlo(arma::mat& X, arma::vec clr_mu, arma::mat clr_sigma, arma::mat& Z){

  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;
  int ds = d;
  int nsim = Z.n_cols;

  arma::vec clr_eigval;
  arma::mat clr_eigvec;
  eig_sym(clr_eigval, clr_eigvec, clr_sigma);
  arma::mat B = clr_eigvec.cols(arma::find(clr_eigval > 1e-5)).t();
  ds = B.n_rows;

  arma::vec B_mu = B * clr_mu;
  arma::mat B_sigma = B * clr_sigma * B.t();
  arma::mat inv_B_sigma = arma::inv_sympd(B_sigma);
  arma::mat Z_ds = Z.head_rows(ds);

  arma::mat Bs_N_mu = arma::zeros(d,n);
  arma::cube Bs_N_sigma = arma::zeros(d,d,n);

  arma::mat clr_E1(D,n);
  arma::cube clr_E2(D,D,n);

  int current_iter = 0;
  arma::vec h(ds);
  arma::vec deriv1(ds);
  arma::mat deriv2(ds, ds);
  arma::vec step(ds);
  for(int k = 0; k<n; k++){
    double eps = 1e-5;
    int max_iter = 100;
    h = B_mu;
    do{
      current_iter++;

      arma::vec Bh = B.t() * h;
      Bh = Bh - Bh.max();

      arma::vec eBh = exp(Bh);

      double  w = arma::accu(eBh);

      arma::vec  wBi = B * eBh;
      deriv1 = -inv_B_sigma  * (h-B_mu) + B * X.col(k) - sum(X.col(k)) * wBi / w;

      arma::mat wBij(ds,ds);
      for(int i=0; i<ds; i++){
        for(int j=0; j<ds; j++){
          wBij(i,j) = arma::accu(B.row(i).t() % B.row(j).t() % eBh);
        }
      }
      deriv2 = -inv_B_sigma - sum(X.col(k)) * ( -(wBi * wBi.t())/ (w*w) + wBij / w);

      step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
      h = h - 0.9 * step;

    }while( norm(step, 2) > eps && current_iter < max_iter);
    arma::vec B_N_mu = h;
    arma::mat B_N_inv_sigma = -deriv2;
    arma::mat B_N_sigma = arma::inv_sympd(B_N_inv_sigma);

    Bs_N_mu.col(k).subvec(0, ds-1) = B_N_mu;
    Bs_N_sigma.slice(k).submat(0, 0, ds-1, ds-1) = B_N_sigma;

    arma::mat B_Z = chol(B_N_sigma).t() * Z_ds;
    B_Z.each_col() += B_N_mu;

    double l_cmult = l_multinomial_const(X.col(k));
    arma::vec l_dens(nsim);
    double l_dens_max = 0;
    double M0 = 0;
    arma::vec M1 = arma::zeros(ds);
    arma::mat M2 = arma::zeros(ds,ds);
    for(int i=0; i<nsim; i++){
      h = B_Z.col(i);
      arma::vec Bh = B.t() * h;
      Bh = Bh - Bh.max();

      arma::vec p = exp(Bh);
      l_dens(i) = l_dnormal_vec(h, B_mu, inv_B_sigma) -
        l_dnormal_vec(h, B_N_mu, B_N_inv_sigma) +
        arma::dot(log(p/accu(p)), X.col(k)) + l_cmult;
      if(l_dens(i) > l_dens_max) l_dens_max = l_dens(i);
    }
    arma::vec dens = exp(l_dens-l_dens_max);
    for(int i=0; i<nsim; i++){
      h = B_Z.col(i);
      M0 += dens(i);
      M1 += h * dens(i);
      M2 += h * h.t() * dens(i);
    }
    clr_E1.col(k) = B.t() * (M1/M0);
    clr_E2.slice(k) = B.t() * (M2/M0) * B;
  }
  List res;
  res["clr_E1"] = clr_E1;
  res["clr_E2"] = clr_E2;
  res["B"] = B;
  res["Bs_N_mu"] = Bs_N_mu;
  res["Bs_N_sigma"] = Bs_N_sigma;
  res["ds"] = ds;
  return(res);
}

//' @export
 // [[Rcpp::export]]
 List c_lrnm_fit_montecarlo(arma::mat& X, arma::mat& Z, double em_eps, int em_max_iter){
   List dm_ini = c_dm_fit(X, 0.001, 1000);
   arma::vec alpha = dm_ini["alpha"];
   arma::mat P(X);
   P.each_col() += alpha;

   int n = X.n_cols;
   int D = X.n_rows;
   int d = D - 1;
   int ds = d;
   int nsim = Z.n_cols;

   arma::mat clr_H = log(P);
   arma::rowvec clr_m = arma::mean(clr_H, 0);
   for(int k=0; k<n; k++) clr_H.col(k) -= clr_m(k);

   arma::vec clr_mu = arma::mean(clr_H, 1);
   arma::mat clr_sigma = arma::cov(clr_H.t());

   arma::mat Bs_N_mu = arma::zeros(d,n);
   arma::cube Bs_N_sigma = arma::zeros(d,d,n);

   arma::mat clr_E1(D,n);
   arma::cube clr_E2(D,D,n);

   arma::vec clr_mu_new;
   arma::mat clr_sigma_new;

   int em_iter = 0;
   double em_current_eps;
   do{
     em_iter++;

     arma::vec clr_eigval;
     arma::mat clr_eigvec;
     eig_sym(clr_eigval, clr_eigvec, clr_sigma);
     arma::mat B = clr_eigvec.cols(arma::find(clr_eigval > 1e-5)).t();
     ds = B.n_rows;

     arma::vec B_mu = B * clr_mu;
     arma::mat B_sigma = B * clr_sigma * B.t();
     arma::mat inv_B_sigma = arma::inv_sympd(B_sigma);
     arma::mat Z_ds = Z.head_rows(ds);


     arma::vec h(ds);
     arma::vec deriv1(ds);
     arma::mat deriv2(ds, ds);
     arma::vec step(ds);

     for(int k = 0; k<n; k++){
       double eps = 1e-5;
       int max_iter = 100;
       h = B * clr_H.col(k);
       int laplace_iter = 0;
       do{
         laplace_iter++;

         arma::vec Bh = B.t() * h;
         Bh = Bh - Bh.max();

         arma::vec eBh = exp(Bh);

         double  w = arma::accu(eBh);

         arma::vec  wBi = B * eBh;
         deriv1 = -inv_B_sigma  * (h-B_mu) + B * X.col(k) - sum(X.col(k)) * wBi / w;

         arma::mat wBij(ds,ds);
         for(int i=0; i<ds; i++){
           for(int j=0; j<ds; j++){
             wBij(i,j) = arma::accu(B.row(i).t() % B.row(j).t() % eBh);
           }
         }
         deriv2 = -inv_B_sigma - sum(X.col(k)) * ( -(wBi * wBi.t())/ (w*w) + wBij / w);

         step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
         h = h - 0.9 * step;

       }while( norm(step, 2) > eps && laplace_iter < max_iter);
       arma::vec B_N_mu = h;
       arma::mat B_N_inv_sigma = -deriv2;
       arma::mat B_N_sigma = arma::inv_sympd(B_N_inv_sigma);

       Bs_N_mu.col(k).subvec(0, ds-1) = B_N_mu;
       Bs_N_sigma.slice(k).submat(0, 0, ds-1, ds-1) = B_N_sigma;

       arma::mat B_Z = chol(B_N_sigma).t() * Z_ds;
       B_Z.each_col() += B_N_mu;

       double l_cmult = l_multinomial_const(X.col(k));
       arma::vec l_dens(nsim);
       double l_dens_max = 0;
       double M0 = 0;
       arma::vec M1 = arma::zeros(ds);
       arma::mat M2 = arma::zeros(ds,ds);
       for(int i=0; i<nsim; i++){
         h = B_Z.col(i);
         arma::vec Bh = B.t() * h;
         Bh = Bh - Bh.max();

         arma::vec p = exp(Bh);
         l_dens(i) = l_dnormal_vec(h, B_mu, inv_B_sigma) -
           l_dnormal_vec(h, B_N_mu, B_N_inv_sigma) +
           arma::dot(log(p/accu(p)), X.col(k)) + l_cmult;
         if(l_dens(i) > l_dens_max) l_dens_max = l_dens(i);
       }
       arma::vec dens = exp(l_dens-l_dens_max);
       for(int i=0; i<nsim; i++){
         h = B_Z.col(i);
         M0 += dens(i);
         M1 += h * dens(i);
         M2 += h * h.t() * dens(i);
       }
       clr_E1.col(k) = B.t() * (M1/M0);
       clr_E2.slice(k) = B.t() * (M2/M0) * B;
     }
     clr_mu_new = mean(clr_E1, 1);
     em_current_eps = norm(clr_mu-clr_mu_new, 1);

     clr_mu = clr_mu_new;
     arma::mat M = mean(clr_E2, 2);
     clr_sigma = M - clr_mu_new * clr_mu_new.t();

   } while (em_iter < em_max_iter & em_current_eps > em_eps);

   List res;
   res["clr_mu"] = clr_mu;
   res["clr_sigma"] = clr_sigma;
   res["ds"] = ds;
   res["clr_E1"] = clr_E1;
   res["em_iter"] = em_iter;
   return(res);
 }
