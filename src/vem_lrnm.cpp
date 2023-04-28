// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include "hermite.h"
#include "dm.h"
#include "lrnm_utils.h"
#include "lrnm_gaussian_approx.h"
#include "coda_base.h"

using namespace Rcpp;


double c_F(arma::vec x, arma::vec m, arma::mat V, arma::vec mu, arma::mat sigma,
           double xi, arma::mat B){

  double K = x.n_elem-1;
  arma::vec mV = m + V.diag()/2;
  arma::mat inv_sigma = arma::pinv(sigma);

  double S1 = arma::dot(x, m);
  double S2 = - arma::accu(x) * ((1/xi) * arma::accu(arma::exp(mV)) - 1 + log(xi));
  double S3 = - 0.5 * log(arma::det(B.t() * sigma * B));
  double S4 = - 0.5 * ((m-mu).t() * inv_sigma * (m-mu)).max();
  double S5 = -0.5 * arma::trace(inv_sigma * V);
  double S6 = 0.5 * (arma::accu(log(V.submat(0,0,K-1,K-1).diag())) + K);
  // Rcpp::Rcout << "S1:" << S1 << std::endl;
  // Rcpp::Rcout << "S2:" << S2 << std::endl;
  // Rcpp::Rcout << "S3:" << S3 << std::endl;
  // Rcpp::Rcout << "S4:" << S4 << std::endl;
  // Rcpp::Rcout << "S5:" << S5 << std::endl;
  // Rcpp::Rcout << "S6:" << S6 << std::endl;
  return(S1 + S2 + S3 + S4 + S5 + S6);
}


// [[Rcpp::export]]
List c_vem_lrnm_fit(arma::mat& X, double em_eps, int em_max_iter){

  List dm_ini = c_dm_fit(X, 0.001, 1000);
  arma::vec alpha = dm_ini["alpha"];
  arma::mat P(X);
  P.each_col() += alpha;

  // Better to condition P on X, C and alpha

  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;


  // int dmax = std::min(d,n-1);
  //
  arma::mat alr_H = log(P.head_rows(d));
  alr_H.each_row() -= log(P.tail_rows(1));

  arma::vec mu = arma::mean(alr_H, 1);
  arma::mat sigma = arma::cov(alr_H.t());

  arma::vec mu_new(d);
  arma::mat sigma_new(d,d);

  arma::mat m = arma::zeros(D,n);
  arma::cube V = arma::zeros(D,D,n);
  for(int k=0;k<n;k++){
    for(int i=0;i<d;i++){
      V(i,i,k) = 1;
    }
  }
  arma::vec xi = arma::zeros(n);
  arma::mat B = arma::eye(D,d);

  arma::mat inv_sigma = arma::inv_sympd(sigma);

  arma::vec mu_ = arma::zeros(D);
  arma::mat sigma_ = arma::zeros(D,D);
  arma::mat inv_sigma_ = arma::zeros(D,D);

  int em_iter = 0;
  double em_current_eps = 1;
  do{
    em_iter++;
    for(int k=0; k<n; k++){
      double x_accu = arma::accu(X.col(k));

      mu_.subvec(0,d-1) = mu;
      sigma_.submat(0,0,d-1,d-1) = sigma;
      inv_sigma_.submat(0,0,d-1,d-1) = arma::inv_sympd(sigma);

      // double logLik = 0;
      // for(int k=0;k<n;k++){
      //   logLik += c_F(X.col(k), m.col(k), V.slice(k), mu_, sigma_, xi(k), B);
      // }
      // Rcpp::Rcout << logLik << std::endl;
      arma::vec mv = m.col(k)+0.5*V.slice(k).diag();
      double mv_max = mv.max();
      xi(k) = arma::accu(arma::exp(mv-mv_max));

      arma::vec m_h = m.col(k);
      arma::vec deriv1(D);
      arma::mat deriv2(D, D);
      arma::vec step = arma::ones(d);
      double eps = 1e-5;
      int max_iter = 400;
      int current_iter=0;
      do{
        current_iter++;
        arma::vec Emv =  arma::exp(m_h + 0.5 * V.slice(k).diag()-mv_max);
        deriv1 = X.col(k) - inv_sigma_ * (m_h-mu_) - x_accu / xi(k) * Emv;
        deriv2 = -inv_sigma_ - x_accu / xi(k) * arma::diagmat(Emv);

        step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
        step(d) = 0;
        m_h = m_h - 0.5 * step;

      }while( norm(step, 2) > eps && current_iter < max_iter);
      if(current_iter >= max_iter) Rcpp::Rcout << "More iterations are needed";
      m.col(k) = m_h;
      for(int i=0;i<d;i++){
        double d1 = 0;
        double d2 = 0;
        double vk = std::sqrt(V(i,i,k));
        double vstep = 1;
        current_iter = 0;
        do{
          current_iter++;
          d1 = 1/vk - vk * inv_sigma_(i,i) - x_accu / xi(k) * exp(m(i,k) + 0.5 * vk*vk-mv_max) * vk;
          d2 = -1/(vk*vk) - x_accu / xi(k) * exp(m(i,k) + 0.5 * vk * vk-mv_max) * (vk*vk +1);
          vstep = d1/d2;
          vk = vk - 0.5 * vstep;
        } while (vstep > eps && current_iter < max_iter);
        if(current_iter >= max_iter) Rcpp::Rcout << "More iterations are needed";
        V(i,i,k) = vk*vk;
      }
    }
    mu_new = mean(m.head_rows(d), 1);
    em_current_eps = norm(mu-mu_new, 1);

    mu = mu_new;
    arma::mat V_mean = mean(V, 2);
    arma::mat EE = arma::zeros(d,d);
    for(int k=0;k<n;k++){
      EE += (m.col(k).head(d)-mu) * (m.col(k).head(d)-mu).t();
    }
    sigma = V_mean.submat(0,0,d-1,d-1) + EE / n;


  } while (em_iter < em_max_iter & em_current_eps > em_eps);

  arma::mat Bclr = (B + 1/D);

  List M;
  M["B"] = B;
  M["alr_h"] = alr_H;
  M["m"] = m;
  M["V"] = V;
  M["xi"] = xi;
  M["clr_H"] = Bclr * m.head_rows(d);
  M["clr_mu"] = Bclr * mu;
  M["clr_sigma"] = Bclr * sigma * Bclr.t();
  M["em_iter"] = em_iter;
  return(M);
}
