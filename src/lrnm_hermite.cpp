// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "dm.h"
#include "hermite.h"
#include "lrnm_utils.h"
#include "lrnm_gaussian_approx.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec c_d_lrnm_hermite_mat(arma::mat X,
                               arma::vec mu_prior, arma::mat sigma_prior,
                               arma::mat Binv, int order){

  unsigned d = X.n_rows - 1;
  unsigned n = X.n_cols;


  arma::mat inv_sigma_prior = arma::inv_sympd(sigma_prior);
  arma::vec dens = arma::zeros(n);
  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  for(int i=0;i<n;i++){
    arma::mat N_posterior = c_lrnm_posterior_approximation_vec(X.col(i), mu_prior, inv_sigma_prior, Binv, 1e-5, 1000);
    arma::vec mu = N_posterior.col(d);
    arma::mat sigma = N_posterior.head_cols(d);
    arma::mat inv_sigma = arma::inv_sympd(sigma);


    arma::vec eigval;
    arma::mat eigvec;

    eig_sym(eigval, eigvec, sigma);
    arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

    unsigned int index[d+1];
    for(unsigned int j = 0; j <= d; j++) index[j] = 0;
    int position = 0, k = 0;

    double l_cmult = l_multinomial_const(X.col(i));

    do{
      double w = 0;
      arma::vec h(d);
      for(unsigned int j = 0; j < d; j++){
        h(j) = uni_hermite(index[j],0);
        w += uni_hermite(index[j],1);
      }
      h = mu + rotation * h;
      arma::vec p = exp(Binv * h);
      dens(i) += exp(w + l_dnormal_vec(h, mu_prior, inv_sigma_prior) -
        l_dnormal_vec(h, mu, inv_sigma) +
        l_multinomial(X.col(i), p/accu(p), l_cmult));
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
  }
  return dens;
}

// [[Rcpp::export]]
List c_lrnm_posterior_moments_hermite(arma::mat& X, arma::vec clr_mu, arma::mat clr_sigma, int order){

  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;
  int ds = d;

  arma::vec clr_eigval;
  arma::mat clr_eigvec;
  eig_sym(clr_eigval, clr_eigvec, clr_sigma);
  arma::mat B = clr_eigvec.cols(arma::find(clr_eigval > 1e-5)).t();
  ds = B.n_rows;

  arma::vec B_mu = B * clr_mu;
  arma::mat B_sigma = B * clr_sigma * B.t();
  arma::mat inv_B_sigma = arma::inv_sympd(B_sigma);

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

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

    double l_cmult = l_multinomial_const(X.col(k));


    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, B_N_sigma);
    arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

    unsigned int index[ds+1];
    for(unsigned int i = 0; i <= ds; i++) index[i] = 0;
    int position = 0, kappa = 0;
    arma::vec l_dens(std::pow(order, ds));
    double l_dens_max = 0;
    do{
      double w = 0;
      for(unsigned int i = 0; i < ds; i++){
        h(i) = uni_hermite(index[i],0);
        w += uni_hermite(index[i],1);
      }
      h = B_N_mu + rotation * h;
      arma::vec Bh = B.t() * h;
      Bh = Bh - Bh.max();

      arma::vec p = exp(Bh);
      l_dens(kappa) = w + l_dnormal_vec(h, B_mu, inv_B_sigma) -
        l_dnormal_vec(h, B_N_mu, B_N_inv_sigma) +
        arma::dot(log(p/accu(p)), X.col(k)) + l_cmult;

      if(l_dens(kappa) > l_dens_max) l_dens_max = l_dens(kappa);


      // Calculate next coordinate
      index[position]++;
      while(index[position] == order){
        index[position] = 0;
        position++;
        index[position]++;
      }
      position = 0;
      kappa++;
    } while (index[ds] == 0);
    arma::vec dens = exp(l_dens-l_dens_max);

    double M0 = 0;
    arma::vec M1 = arma::zeros(ds);
    arma::mat M2 = arma::zeros(ds,ds);

    for(unsigned int i = 0; i <= ds; i++) index[i] = 0;
    position = 0;
    kappa = 0;
    do{
      for(unsigned int i = 0; i < ds; i++){
        h(i) = uni_hermite(index[i],0);
      }
      h = B_N_mu + rotation * h;

      M0 += dens(kappa);
      M1 += h * dens(kappa);
      M2 += h * h.t() * dens(kappa);

      // Calculate next coordinate
      index[position]++;
      while(index[position] == order){
        index[position] = 0;
        position++;
        index[position]++;
      }
      position = 0;
      kappa++;
    } while (index[ds] == 0);

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


// [[Rcpp::export]]
List c_lrnm_fit_hermite(arma::mat& X, int order, double em_eps, int em_max_iter){
  List dm_ini = c_dm_fit(X, 0.001, 1000);
  arma::vec alpha = dm_ini["alpha"];
  arma::mat P(X);
  P.each_col() += alpha;

  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;
  int ds = d;

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::mat clr_H = log(P);
  arma::rowvec clr_m = arma::mean(clr_H, 0);
  for(int k=0; k<n; k++) clr_H.col(k) -= clr_m(k);

  arma::vec clr_mu = arma::mean(clr_H, 1);
  arma::mat clr_sigma = arma::cov(clr_H.t());

  arma::mat B;
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
    B = clr_eigvec.cols(arma::find(clr_eigval > 1e-5)).t();
    ds = B.n_rows;

    arma::vec B_mu = B * clr_mu;
    arma::mat B_sigma = B * clr_sigma * B.t();
    arma::mat inv_B_sigma = arma::inv_sympd(B_sigma);

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


      double l_cmult = l_multinomial_const(X.col(k));

      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, B_N_sigma);
      arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

      unsigned int index[ds+1];
      for(unsigned int i = 0; i <= ds; i++) index[i] = 0;
      int position = 0, kappa = 0;
      arma::vec l_dens(std::pow(order, ds));
      double l_dens_max = 0;
      do{
        double w = 0;
        for(unsigned int i = 0; i < ds; i++){
          h(i) = uni_hermite(index[i],0);
          w += uni_hermite(index[i],1);
        }
        h = B_N_mu + rotation * h;
        arma::vec Bh = B.t() * h;
        Bh = Bh - Bh.max();

        arma::vec p = exp(Bh);
        l_dens(kappa) = w + l_dnormal_vec(h, B_mu, inv_B_sigma) -
          l_dnormal_vec(h, B_N_mu, B_N_inv_sigma) +
          arma::dot(log(p/accu(p)), X.col(k)) + l_cmult;

        if(l_dens(kappa) > l_dens_max) l_dens_max = l_dens(kappa);


        // Calculate next coordinate
        index[position]++;
        while(index[position] == order){
          index[position] = 0;
          position++;
          index[position]++;
        }
        position = 0;
        kappa++;
      } while (index[ds] == 0);
      arma::vec dens = exp(l_dens-l_dens_max);

      double M0 = 0;
      arma::vec M1 = arma::zeros(ds);
      arma::mat M2 = arma::zeros(ds,ds);

      for(unsigned int i = 0; i <= ds; i++) index[i] = 0;
      position = 0;
      kappa = 0;
      do{
        for(unsigned int i = 0; i < ds; i++){
          h(i) = uni_hermite(index[i],0);
        }
        h = B_N_mu + rotation * h;

        M0 += dens(kappa);
        M1 += h * dens(kappa);
        M2 += h * h.t() * dens(kappa);

        // Calculate next coordinate
        index[position]++;
        while(index[position] == order){
          index[position] = 0;
          position++;
          index[position]++;
        }
        position = 0;
        kappa++;
      } while (index[ds] == 0);

      clr_E1.col(k) = B.t() * (M1/M0);
      clr_E2.slice(k) = B.t() * (M2/M0) * B;


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
  res["Bs"] = B;
  res["Bs_N_mu"] = Bs_N_mu;
  res["Bs_N_sigma"] = Bs_N_sigma;
  return(res);
}


// [[Rcpp::export]]
List c_lrnm_fit_fixed_hermite(arma::mat& X, int order, double em_eps, int em_max_iter){
  List dm_ini = c_dm_fit(X, 0.001, 1000);
  arma::vec alpha = dm_ini["alpha"];
  arma::mat P(X);
  P.each_col() += alpha;

  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;

  arma::mat clr_H = log(P);
  arma::rowvec clr_m = arma::mean(clr_H, 0);
  for(int k=0; k<n; k++) clr_H.col(k) -= clr_m(k);

  arma::vec clr_mu = arma::mean(clr_H, 1);
  arma::mat clr_sigma = arma::cov(clr_H.t());

  if(em_max_iter == 0){
    List res;
    res["clr_H"] = clr_H;
    res["clr_mu"] = clr_mu;
    res["clr_sigma"] = clr_sigma;
    return(res);
  }

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::vec clr_eigval;
  arma::mat clr_eigvec;
  eig_sym(clr_eigval, clr_eigvec, clr_sigma);

  arma::mat B = clr_eigvec.cols(arma::find(clr_eigval > 1e-8)).t();

  int ds = B.n_rows;

  arma::vec B_mu(ds);
  arma::mat B_sigma(ds,ds);
  arma::mat B_inv_sigma(ds,ds);

  arma::vec B_N_mu(ds);
  arma::mat B_N_sigma(ds,ds);
  arma::mat B_inv_N_sigma(ds,ds);

  arma::mat clr_E1(clr_H);
  arma::cube clr_E2(D,D,n);

  arma::vec clr_mu_new(D);
  arma::mat clr_sigma_new(D,D);

  arma::mat Ieps = 1e-8 * arma::eye(ds,ds);
  arma::vec l_cmult(n), Xsum(n);
  for(int k = 0; k<n; k++){
    l_cmult(k) = l_multinomial_const(X.col(k));
    Xsum(k) = sum(X.col(k));
  }

  int em_iter = 0;
  double em_current_eps;
  do{
    Rcpp::checkUserInterrupt();

    em_iter++;

    B_mu = B * clr_mu;
    B_sigma = Ieps + B * clr_sigma * B.t();
    B_inv_sigma = arma::inv_sympd(B_sigma);

    arma::vec h(ds);
    arma::vec deriv1(ds);
    arma::mat deriv2(ds, ds);
    arma::vec step(ds);
    double n_step;
    arma::mat H = B * clr_E1;
    for(int k = 0; k<n; k++){
      double eps = 1e-6;
      int max_iter = 100;
      h = H.col(k);
      double nMax = norm(h - B_mu);
      int laplace_iter = 0;
      do{
        laplace_iter++;

        arma::vec Bh = B.t() * h;
        Bh = Bh - Bh.max();

        arma::vec eBh = exp(Bh);

        double  w = arma::accu(eBh);

        arma::vec  wBi = B * eBh;
        deriv1 = -B_inv_sigma  * (h-B_mu) + B * X.col(k) - Xsum(k) * wBi / w;

        arma::mat wBij(ds,ds);
        for(int i=0; i<ds; i++){
          for(int j=0; j<ds; j++){
            wBij(i,j) = arma::accu(B.row(i).t() % B.row(j).t() % eBh);
          }
        }
        deriv2 = -B_inv_sigma - Xsum(k) * ( -(wBi * wBi.t())/ (w*w) + wBij / w);

        step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
        n_step = norm(step, 2);
        if(n_step < 0.5*nMax){
          h = h - step;
        }else{
          h = h - 0.5*nMax/n_step * step;
        }

      }while( n_step > eps && laplace_iter < max_iter );

      B_N_mu = h;
      B_inv_N_sigma = -deriv2;
      B_N_sigma = arma::inv_sympd(B_inv_N_sigma);


      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, B_N_sigma);
      arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

      unsigned int index[ds+1];
      for(unsigned int i = 0; i <= ds; i++) index[i] = 0;
      int position = 0, kappa = 0;
      arma::vec l_dens(std::pow(order, ds));
      do{
        double w = 0;
        for(unsigned int i = 0; i < ds; i++){
          h(i) = uni_hermite(index[i],0);
          w += uni_hermite(index[i],1);
        }
        h = B_N_mu + rotation * h;
        arma::vec Bh = B.t() * h;
        Bh = Bh - Bh.max();

        arma::vec p = exp(Bh);
        l_dens(kappa) = w + l_dnormal_vec(h, B_mu, B_inv_sigma) -
          l_dnormal_vec(h, B_N_mu, B_inv_N_sigma) +
          arma::dot(Bh, X.col(k)) - accu(log(accu(p)) * X.col(k)) + l_cmult(k);
          // arma::dot(log(p/accu(p)), X.col(k)) + l_cmult(k);

        // Calculate next coordinate
        index[position]++;
        while(index[position] == order){
          index[position] = 0;
          position++;
          index[position]++;
        }
        position = 0;
        kappa++;
      } while (index[ds] == 0);
      arma::vec dens = exp(l_dens-l_dens.max());

      double M0 = 0;
      arma::vec M1 = arma::zeros(ds);
      arma::mat M2 = arma::zeros(ds,ds);

      for(unsigned int i = 0; i <= ds; i++) index[i] = 0;
      position = 0;
      kappa = 0;
      do{
        for(unsigned int i = 0; i < ds; i++){
          h(i) = uni_hermite(index[i],0);
        }
        h = B_N_mu + rotation * h;

        M0 += dens(kappa);
        M1 += h * dens(kappa);
        M2 += h * h.t() * dens(kappa);

        // Calculate next coordinate
        index[position]++;
        while(index[position] == order){
          index[position] = 0;
          position++;
          index[position]++;
        }
        position = 0;
        kappa++;
      } while (index[ds] == 0);

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
  res["clr_H"] = clr_E1;
  res["clr_mu"] = clr_mu;
  res["clr_sigma"] = clr_sigma;
  res["B"] = B;
  res["d"] = ds;
  res["em_iter"] = em_iter;
  res["em_final_eps"] = em_current_eps;
  return(res);
}
