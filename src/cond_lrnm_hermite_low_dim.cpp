// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include "dm.h"
#include "coda_base.h"
#include "lrnm_utils.h"
#include "hermite.h"

using namespace Rcpp;

// [[Rcpp::export]]
List c_cond_lrnm_one_dim_fit_hermite(arma::mat& X, arma::mat& C,
                                     int order,
                                     double em_eps, int em_max_iter){
  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;

  arma::vec d0 = arma::vec(n);
  arma::umat I0 = arma::umat(D, n);
  for(int k=0; k<n; k++){
    d0(k) = arma::accu(C.col(k) == 0);
    I0.col(k).head(d0(k)) = arma::find(C.col(k) == 0);
  }

  arma::vec d1 = d-d0;
  arma::vec D1 = 1+d1;
  arma::umat I1 = arma::umat(D, n);
  for(int k=0; k<n; k++){
    I1.col(k).head(D1(k)) = arma::find(C.col(k) > 0);
  }

  List lrnm_res = c_dm_fit(X, 0.001, 1000);
  arma::vec alpha = lrnm_res["alpha"];
  arma::mat P(X);
  for(unsigned int k=0;k<n;k++){
    arma::uvec uk(1);
    uk[0] = k;
    arma::vec a_p = X.col(k) + alpha;
    arma::vec a_p_1 = a_p / arma::accu(a_p);
    double res = arma::accu(a_p_1(I1.col(k).head(D1(k))));
    arma::vec a_p_0 = a_p(I0.col(k).head(d0(k)));
    a_p_0 = a_p_0 / arma::accu(a_p_0) * (1-res);
    P.col(k) = P.col(k) / arma::accu(P.col(k)) * res;
    P(I0.col(k).head(d0(k)), uk) = a_p_0;
  }


  arma::mat clr_H = log(P);
  arma::rowvec clr_m = arma::mean(clr_H, 0);
  for(int k=0; k<n; k++) clr_H.col(k) -= clr_m(k);

  arma::vec clr_mu = arma::mean(clr_H, 1);
  arma::mat clr_sigma = arma::cov(clr_H.t());

  arma::cube B1 = arma::zeros(D,d,n);
  arma::vec d1s = arma::zeros(n);
  arma::vec d1_unique = unique(d1);
  arma::cube ILR1 = arma::zeros(D,d,d1_unique.n_elem);

  arma::vec i1 = arma::zeros(D+1);
  for(int i=0;i<d1_unique.n_elem;i++){
    i1(d1_unique(i)) = i;
    if(d1_unique(i) > 0) ILR1.slice(i).submat(0,0,d1_unique(i),d1_unique(i)-1) = ilr_basis_default(d1_unique(i)+1);
  }
  for(int k=0;k<n;k++){
    if(d1(k) > 0){
      arma::uvec rows = I1.col(k).head(D1(k));
      arma::mat ilr_m = ILR1.slice(i1(d1(k))).head_rows(D1(k));
      B1.slice(k).rows(rows) = ilr_m;

      arma::mat B1s = B1.slice(k).head_cols(d1(k));
      arma::mat B1_sigma = B1s.t() * clr_sigma * B1s;

      arma::vec B1_eigval;
      arma::mat B1_eigvec;
      eig_sym(B1_eigval, B1_eigvec, B1_sigma);
      arma::mat B1_ = B1_eigvec.cols(arma::find(B1_eigval > 1e-3));
      d1s(k) = B1_.n_cols;
      B1.slice(k).head_cols(d1s(k)) = B1.slice(k).head_cols(d1(k)) * B1_;
    }
  }


  arma::cube B0 = arma::zeros(D,d,n);
  arma::vec d0s = arma::zeros(n);
  arma::vec d0_unique = unique(d0);
  arma::cube ILR0 = arma::zeros(D,d,d0_unique.n_elem);

  arma::vec i0 = arma::zeros(D+1);
  for(int i=0;i<d0_unique.n_elem;i++){
    i0(d0_unique(i)) = i;
    if(d0_unique(i)>1){
      ILR0.slice(i).submat(0,0,d0_unique(i)-1,d0_unique(i)-2) = ilr_basis_default(d0_unique(i));
    }
  }
  for(int k=0;k<n;k++){

    if(d0(k) > 0){
      if(d0(k) > 1){
        arma::uvec rows = I0.col(k).head(d0(k));
        arma::mat ilr_m = ILR0.slice(i0(d0(k))).head_rows(d0(k));
        B0.slice(k).rows(rows) = ilr_m;
      }
      double v0 = +1/d0(k) * std::sqrt(d0(k) * (D1(k))/D);
      for(int i=0; i<d0(k);i++){
        B0.slice(k)(I0(i,k),d0(k)-1) = v0;
      }
      double v1 = -1/(D1(k)) * std::sqrt(d0(k) * (D1(k))/D);
      for(int i=0; i<D1(k);i++){
        B0.slice(k)(I1(i,k),d0(k)-1) = v1;
      }
      arma::mat B0s = B0.slice(k).head_cols(d0(k));
      arma::mat B0_sigma = B0s.t() * clr_sigma * B0s;
      // Basis reduction
      if(d1(k)>0){
        arma::mat B1s = B1.slice(k).head_cols(d1s(k));
        arma::vec h1 = B1s.t() * clr_H.col(k);
        arma::mat clr_sigma_B1s = clr_sigma * B1s;
        arma::mat B1_0 = B0s.t() * clr_sigma_B1s * arma::inv_sympd(B1s.t() * clr_sigma_B1s);
        // B1_0.print("B1_0");
        // B0_sigma.print("B0_sigma");
        B0_sigma = B0_sigma - B1_0 * clr_sigma_B1s.t() * B0s;
        // B0_sigma.print("B0_sigma");
      }
      arma::vec B0_eigval;
      arma::mat B0_eigvec;
      eig_sym(B0_eigval, B0_eigvec, B0_sigma);
      arma::mat B0_ = B0_eigvec.cols(arma::find(B0_eigval > 1e-3));
      d0s(k) = B0_.n_cols;
      B0.slice(k).head_cols(d0s(k)) = B0.slice(k).head_cols(d0(k)) * B0_;
    }
  }

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::mat B0v = arma::zeros(d,n);
  arma::mat V_N_mu(1,n);
  arma::cube V_N_sigma(1,1,n);

  arma::mat EV1(1,n);
  arma::cube EV2(1,1,n);

  arma::mat clr_E1(clr_H);
  arma::cube clr_E2(D,D,n);
  for(int k=0;k<n;k++) clr_E2.slice(k) = clr_E1.col(k) * clr_E1.col(k).t();

  ///
  int em_iter = 0;
  double em_current_eps;
  arma::vec clr_mu_new;
  do{
    em_iter++;
    for(int k=0;k<n;k++){
      if(d0s(k)>0){
        arma::mat B0s = B0.slice(k).head_cols(d0s(k));
        arma::mat B1s = B1.slice(k).head_cols(d1s(k));
        arma::vec clr_h = clr_E1.col(k);
        arma::vec h1 = B1s.t() * clr_h;
        arma::vec h0 = B0s.t() * clr_h;

        // (B1s.t() * clr_sigma * B1s).print("B1s.t() * clr_sigma * B1s");
        arma::mat inv_B1_sigma = arma::inv_sympd(B1s.t() * clr_sigma * B1s);
        // inv_B1_sigma.print("inv_B1");
        arma::mat B1_0 = B1s.t() * B0s;

        arma::vec B0_mu = B0s.t() * clr_mu + B1_0.t() * inv_B1_sigma * (h1 - B1s.t() * clr_mu);
        arma::mat B0_sigma = B0s.t() * clr_sigma * B0s - B1_0.t() * inv_B1_sigma * B1_0;
        arma::mat inv_B0_sigma = arma::inv_sympd(B0_sigma);


        arma::vec V = arma::normalise(h0 - B0_mu);

        B0v.col(k).head(d0s(k)) = V;

        int current_iter = 0;
        int dh0 = 1;
        arma::vec deriv1(1), deriv2(1), step = arma::ones(1);
        arma::vec eps = {1e-5};
        int max_iter = 100;

        arma::vec v = V.t() * h0; // {0}; //
        arma::vec h0_(d0s(k));
        arma::vec h(D);
        do{

          current_iter++;
          h0_ = B0_mu + V * v;
          if(d1s(k)>0){
            h = B0s * h0_ + B1s * h1;
          }else{
            h = B0s * h0_;
          }
          arma::vec eBh = exp(h);
          double  w = arma::accu(eBh);

          arma::vec  wBi = B0s.t() * eBh;
          deriv1 = V.t() * (-inv_B0_sigma  * (h0_-B0_mu) + B0s.t() * X.col(k) - sum(X.col(k)) * wBi / w);

          arma::mat wBij(d0s(k), d0s(k));
          for(int i=0;i<d0s(k);i++){
            for(int j=0;j<d0s(k);j++){
              wBij(i,j) = arma::accu(B0s.col(i) % B0s.col(j) % eBh);
            }
          }

          deriv2 = V.t() * (-inv_B0_sigma - sum(X.col(k)) * ( -(wBi * wBi.t())/ (w*w) + wBij / w)) * V;
          step = deriv1/deriv2;
          v = v - 0.9 * step;

        }while( std::abs(step(0)) > eps(0) && current_iter < max_iter);
        arma::vec v_n_mu = v;
        arma::mat v_n_sigma = -1/deriv2;
        arma::mat inv_v_n_sigma = -deriv2;

        V_N_mu.col(k) = v;
        V_N_sigma.slice(k) = -1/deriv2;

        double l_cmult = l_multinomial_const(X.col(k));

        int dh0_red = 1;
        unsigned int index[] = {0,0};
        int position = 0, kappa = 0;
        arma::vec l_dens(order);

        double w;
        do{
          v = v_n_mu + sqrt(v_n_sigma(0,0)) * uni_hermite(index[0],0);
          w = uni_hermite(index[0],dh0_red);


          h0_ = B0_mu + V * v;
          if(d1s(k)>0){
            h = B0s * h0_ + B1s * h1;
          }else{
            h = B0s * h0_;
          }
          h = h - h.max();
          arma::vec p = exp(h);
          l_dens(kappa) = w + l_dnormal_vec(h0_, B0_mu, inv_B0_sigma) -
            l_dnormal_vec(v, v_n_mu, inv_v_n_sigma) +
            arma::dot(log(p/accu(p)), X.col(k)) + l_cmult;

          // Calculate next coordinate
          index[position]++;
          while(index[position] == order){
            index[position] = 0;
            position++;
            index[position]++;
          }
          position = 0;
          kappa++;
        } while (index[(unsigned)dh0_red] == 0);
        arma::vec dens = exp(l_dens-l_dens.max());

        double M0 = 0;
        arma::vec M1 = arma::zeros(dh0_red);
        arma::mat M2 = arma::zeros(dh0_red,dh0_red);
        for(unsigned int i = 0; i <= dh0_red; i++) index[i] = 0;
        position = 0;
        kappa = 0;
        do{
          v = v_n_mu + sqrt(v_n_sigma(0,0)) * uni_hermite(index[0],0);

          M0 += dens(kappa);
          M1 += v * dens(kappa);
          M2 += v.t() * v * dens(kappa);

          // Calculate next coordinate
          index[position]++;
          while(index[position] == order){
            index[position] = 0;
            position++;
            index[position]++;
          }
          position = 0;
          kappa++;
        } while (index[(unsigned)dh0_red] == 0);
        EV1.col(k) = M1/M0;
        EV2.slice(k) = M2/M0;

        // cbind(B0,B1) %*% rbind(
        //     cbind(V %*% cbind(I2$EV2[,,1]) %*% t(V), V %*% cbind(I2$EV1[,1]) %*% t(h1)),
        //     cbind(h1 %*% cbind(I2$EV1[,1]) %*% t(V), h1 %*% t(h1))) %*% t(cbind(B0,B1))


        arma::vec Eh0 = B0_mu + V * (M1/M0);
        arma::mat Eh0h0 = B0_mu * B0_mu.t() + B0_mu * (V * M1/M0).t() + V * (M1/M0) * B0_mu.t() + V * (M2/M0) * V.t();

        if(d1s(k)>0){
          arma::mat B = arma::join_rows(B0s,B1s);
          clr_E1.col(k) = B * arma::join_cols(Eh0, h1);
          clr_E2.slice(k) = B * arma::join_cols(
            arma::join_rows(Eh0h0, Eh0 * h1.t()),
            arma::join_rows(h1 * Eh0.t(), h1 * h1.t())) * B.t();
        }else{
          clr_E1.col(k) = B0s * Eh0;
          clr_E2.slice(k) = B0s * Eh0h0 * B0s.t();
        }
      }
    }
    clr_mu_new = mean(clr_E1, 1);
    em_current_eps = norm(clr_mu-clr_mu_new, 1);

    clr_mu = clr_mu_new;
    arma::mat M = mean(clr_E2, 2);
    clr_sigma = M - clr_mu_new * clr_mu_new.t();
  }while(em_iter < em_max_iter & em_current_eps > em_eps);

  lrnm_res["d0"] = d0;
  lrnm_res["d1"] = d1;
  lrnm_res["d0s"] = d0s;
  lrnm_res["d1s"] = d1s;
  lrnm_res["I0"] = I0;
  lrnm_res["I1"] = I1;
  lrnm_res["B0"] = B0;
  lrnm_res["B1"] = B1;
  lrnm_res["P"] = P;
  lrnm_res["X"] = X;
  lrnm_res["clr_H"] = clr_E1;
  lrnm_res["clr_mu"] = clr_mu;
  lrnm_res["clr_sigma"] = clr_sigma;
  lrnm_res["B0v"] = B0v;

  lrnm_res["V_N_mu"] = V_N_mu;
  lrnm_res["V_N_sigma"] = V_N_sigma;
  lrnm_res["clr_E1"] = clr_E1;
  lrnm_res["em_iter"] = em_iter;
  return(lrnm_res);
}

// [[Rcpp::export]]
List c_cond_lrnm_full_laplace_one_dim_fit_hermite(arma::mat& X, arma::mat& C,
                                                  int order,
                                                  double em_eps, int em_max_iter){
  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;

  arma::vec d0 = arma::vec(n);
  arma::umat I0 = arma::umat(D, n);
  for(int k=0; k<n; k++){
    d0(k) = arma::accu(C.col(k) == 0);
    I0.col(k).head(d0(k)) = arma::find(C.col(k) == 0);
  }

  arma::vec d1 = d-d0;
  arma::vec D1 = 1+d1;
  arma::umat I1 = arma::umat(D, n);
  for(int k=0; k<n; k++){
    I1.col(k).head(D1(k)) = arma::find(C.col(k) > 0);
  }

  List lrnm_res = c_dm_fit(X, 0.001, 1000);
  arma::vec alpha = lrnm_res["alpha"];
  arma::mat P(X);
  for(unsigned int k=0;k<n;k++){
    arma::uvec uk(1);
    uk[0] = k;
    arma::vec a_p = X.col(k) + alpha;
    arma::vec a_p_1 = a_p / arma::accu(a_p);
    double res = arma::accu(a_p_1(I1.col(k).head(D1(k))));
    arma::vec a_p_0 = a_p(I0.col(k).head(d0(k)));
    a_p_0 = a_p_0 / arma::accu(a_p_0) * (1-res);
    P.col(k) = P.col(k) / arma::accu(P.col(k)) * res;
    P(I0.col(k).head(d0(k)), uk) = a_p_0;
  }


  arma::mat clr_H = log(P);
  arma::rowvec clr_m = arma::mean(clr_H, 0);
  for(int k=0; k<n; k++) clr_H.col(k) -= clr_m(k);

  arma::vec clr_mu = arma::mean(clr_H, 1);
  arma::mat clr_sigma = arma::cov(clr_H.t());

  if(em_max_iter==0){
    List res;
    res["clr_H"] = clr_H;
    res["clr_mu"] = clr_mu;
    res["clr_sigma"] = clr_sigma;
    return(res);
  }

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::cube B = arma::zeros(d,D,n);
  for(int k=0; k<n; k++){
    if(d0(k) > 0){
      if(d0(k) > 1){
        arma::uvec i0 = I0.col(k).head(d0(k));
        for(unsigned int i = 0; i < d0(k) - 1; i++){
          unsigned int I1 = i + 1;
          unsigned int I2 = i + 2;
          double l = 1/std::sqrt((double)(I1*I2));
          double r = - std::sqrt((double)I1/I2);
          for(unsigned int j = 0; j < I1; j++){
            B(i,i0(j),k) = l;
          }
          B(i,i0(I1),k) = r;
        }
      }
      for(unsigned int j = 0; j < D; j++){
        if(C(j,k) == 0){
          B(d0(k) - 1, j, k) = +1/d0(k) * std::sqrt(d0(k) * (D1(k))/D);
        }else{
          B(d0(k) - 1, j, k) = -1/(D1(k)) * std::sqrt(d0(k) * (D1(k))/D);
        }
      }
      if(D1(k) > 1){
        arma::uvec i1 = I1.col(k).head(D1(k));
        for(unsigned int i = 0; i < D1(k) - 1; i++){
          unsigned int I1 = i + 1;
          unsigned int I2 = i + 2;
          double l = 1/std::sqrt((double)(I1*I2));
          double r = - std::sqrt((double)I1/I2);
          for(unsigned int j = 0; j < I1; j++){
            B(d0(k)+i,i1(j),k) = l;
          }
          B(d0(k)+i,i1(I1),k) = r;
        }
      }
    }else{
      B.slice(k) = ilr_basis_default(D1(k)).t();
    }
  }


  arma::cube B0 = arma::zeros(d,D,n);
  arma::vec d0s = arma::zeros(n);

  arma::cube B1 = arma::zeros(d,D,n);
  arma::vec d1s = arma::zeros(n);
  for(int k = 0; k<n; k++){
    if(d0(k) > 0){
      arma::mat B0_ = B.slice(k).head_rows(d0(k));

      arma::vec clr_eigval;
      arma::mat clr_eigvec;
      eig_sym(clr_eigval, clr_eigvec, B0_ * clr_sigma * B0_.t());

      arma::mat B0s_ = clr_eigvec.cols(arma::find(clr_eigval > 1e-8)).t();
      d0s(k) = B0s_.n_rows;

      B0.slice(k).head_rows(d0s(k)) = B0s_ * B0_;
    }
    if(d1(k) > 0){
      arma::mat B1_ = B.slice(k).tail_rows(d1(k));

      arma::vec clr_eigval;
      arma::mat clr_eigvec;
      eig_sym(clr_eigval, clr_eigvec, B1_ * clr_sigma * B1_.t());

      arma::mat B1s_ = clr_eigvec.cols(arma::find(clr_eigval > 1e-8)).t();
      d1s(k) = B1s_.n_rows;

      B1.slice(k).head_rows(d1s(k)) = B1s_ * B1_;
    }
  }

  arma::mat B0s_N_mu = arma::zeros(d,n);
  arma::cube B0s_N_sigma = arma::zeros(d,d,n);

  // arma::cube B1 = arma::zeros(D,d,n);
  // arma::vec d1s = arma::zeros(n);
  // arma::vec d1_unique = unique(d1);
  // arma::cube ILR1 = arma::zeros(D,d,d1_unique.n_elem);
  //
  // arma::vec i1 = arma::zeros(D+1);
  // for(int i=0;i<d1_unique.n_elem;i++){
  //   i1(d1_unique(i)) = i;
  //   if(d1_unique(i) > 0) ILR1.slice(i).submat(0,0,d1_unique(i),d1_unique(i)-1) = ilr_basis_default(d1_unique(i)+1);
  // }
  // for(int k=0;k<n;k++){
  //   if(d1(k) > 0){
  //     arma::uvec rows = I1.col(k).head(D1(k));
  //     arma::mat ilr_m = ILR1.slice(i1(d1(k))).head_rows(D1(k));
  //     B1.slice(k).rows(rows) = ilr_m;
  //
  //     arma::mat B1s = B1.slice(k).head_cols(d1(k));
  //     arma::mat B1_sigma = B1s.t() * clr_sigma * B1s;
  //
  //     arma::vec B1_eigval;
  //     arma::mat B1_eigvec;
  //     eig_sym(B1_eigval, B1_eigvec, B1_sigma);
  //     arma::mat B1_ = B1_eigvec.cols(arma::find(B1_eigval > 1e-3));
  //     d1s(k) = B1_.n_cols;
  //     B1.slice(k).head_cols(d1s(k)) = B1.slice(k).head_cols(d1(k)) * B1_;
  //   }
  // }
  //
  //
  // arma::cube B0 = arma::zeros(D,d,n);
  // arma::vec d0s = arma::zeros(n);
  // arma::vec d0_unique = unique(d0);
  // arma::cube ILR0 = arma::zeros(D,d,d0_unique.n_elem);
  //
  // arma::vec i0 = arma::zeros(D+1);
  // for(int i=0;i<d0_unique.n_elem;i++){
  //   i0(d0_unique(i)) = i;
  //   if(d0_unique(i)>1){
  //     ILR0.slice(i).submat(0,0,d0_unique(i)-1,d0_unique(i)-2) = ilr_basis_default(d0_unique(i));
  //   }
  // }
  // for(int k=0;k<n;k++){
  //
  //   if(d0(k) > 0){
  //     if(d0(k) > 1){
  //       arma::uvec rows = I0.col(k).head(d0(k));
  //       arma::mat ilr_m = ILR0.slice(i0(d0(k))).head_rows(d0(k));
  //       B0.slice(k).rows(rows) = ilr_m;
  //     }
  //     double v0 = +1/d0(k) * std::sqrt(d0(k) * (D1(k))/D);
  //     for(int i=0; i<d0(k);i++){
  //       B0.slice(k)(I0(i,k),d0(k)-1) = v0;
  //     }
  //     double v1 = -1/(D1(k)) * std::sqrt(d0(k) * (D1(k))/D);
  //     for(int i=0; i<D1(k);i++){
  //       B0.slice(k)(I1(i,k),d0(k)-1) = v1;
  //     }
  //     arma::mat B0s = B0.slice(k).head_cols(d0(k));
  //     arma::mat B0_sigma = B0s.t() * clr_sigma * B0s;
  //     // Basis reduction
  //     if(d1(k)>0){
  //       arma::mat B1s = B1.slice(k).head_cols(d1s(k));
  //       arma::vec h1 = B1s.t() * clr_H.col(k);
  //       arma::mat clr_sigma_B1s = clr_sigma * B1s;
  //       arma::mat B1_0 = B0s.t() * clr_sigma_B1s * arma::inv_sympd(B1s.t() * clr_sigma_B1s);
  //       // B1_0.print("B1_0");
  //       // B0_sigma.print("B0_sigma");
  //       B0_sigma = B0_sigma - B1_0 * clr_sigma_B1s.t() * B0s;
  //       // B0_sigma.print("B0_sigma");
  //     }
  //     arma::vec B0_eigval;
  //     arma::mat B0_eigvec;
  //     eig_sym(B0_eigval, B0_eigvec, B0_sigma);
  //     arma::mat B0_ = B0_eigvec.cols(arma::find(B0_eigval > 1e-3));
  //     d0s(k) = B0_.n_cols;
  //     B0.slice(k).head_cols(d0s(k)) = B0.slice(k).head_cols(d0(k)) * B0_;
  //   }
  // }



  arma::mat B0v = arma::zeros(D,n);
  // arma::mat V_N_mu(1,n);
  // arma::cube V_N_sigma(1,1,n);
  //
  // arma::mat EV1(1,n);
  // arma::cube EV2(1,1,n);

  arma::mat clr_E1(clr_H);
  arma::cube clr_E2(D,D,n);
  for(int k=0;k<n;k++) clr_E2.slice(k) = clr_E1.col(k) * clr_E1.col(k).t();

  ///
  int em_iter = 0;
  double em_current_eps;
  arma::vec clr_mu_new;
  do{
    Rcpp::checkUserInterrupt();

    em_iter++;

    for(int k=0;k<n;k++){
      int dh0 = d0s(k);
      int dh1 = d1s(k);

      if(dh0>0){
        arma::mat B0_ = B0.slice(k).head_rows(dh0);
        arma::mat B1_ = B1.slice(k).head_rows(dh1);
        arma::vec B0_mu = B0_ * clr_mu;
        arma::mat B0_sigma = B0_ * clr_sigma * B0_.t();

        if(dh1>0){
          arma::vec h1_mu = B1_ * (clr_E1.col(k)-clr_mu);
          arma::mat B1_sigma = B1_ * clr_sigma * B1_.t();

          arma::mat B1_0 = B0_ * clr_sigma * B1_.t();
          arma::mat B1_0_inv_B1_sigma =  B1_0 * arma::inv_sympd(B1_sigma);

          B0_mu = B0_mu + B1_0_inv_B1_sigma * (h1_mu);
          B0_sigma = B0_sigma - B1_0_inv_B1_sigma * B1_0.t();
        }
        // Finding a Gaussian approximation
        B0_sigma.diag() += 1e-8;
        arma::mat inv_B0_sigma = arma::inv_sympd(B0_sigma);

        arma::vec deriv1(dh0);
        arma::mat deriv2(dh0, dh0);
        arma::vec step = arma::ones(dh0);
        double n_step;
        double eps = 1e-6;
        int max_iter = 100;

        arma::mat Binv = B.slice(k).t();
        arma::vec h0 = B0_ * clr_E1.col(k);
        arma::vec clr_h1 = B1_.t() * B1_ * clr_E1.col(k);
        double nMax = norm(h0 - B0_mu);
        int laplace_iter = 0;
        do{
          laplace_iter++;
          arma::vec clr_h = B0_.t() * h0;
          if(d1(k)>0){
            clr_h += clr_h1;
          }
          arma::vec Bh = clr_h - clr_h.max();
          arma::vec eBh = exp(Bh);

          double  w = arma::accu(eBh);

          arma::vec  wBi = B0_ * eBh;
          deriv1 = -inv_B0_sigma  * (h0-B0_mu) + B0_ * X.col(k) - sum(X.col(k)) * wBi / w;
          arma::mat wBij(dh0,dh0);
          for(int i=0; i<dh0; i++){
            for(int j=0; j<dh0; j++){
              wBij(i,j) = arma::accu(B0_.row(i) % B0_.row(j) % eBh.t());
            }
          }
          deriv2 = -inv_B0_sigma - sum(X.col(k)) * ( -(wBi * wBi.t())/ (w*w) + wBij / w);

          step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
          n_step = norm(step, 2);
          if(n_step < 0.5*nMax){
            h0 = h0 - step;
          }else{
            h0 = h0 - 0.5*nMax/n_step * step;
          }
        }while( n_step > eps && laplace_iter < max_iter );

        arma::vec B0_N_mu = h0;
        arma::mat B0_N_inv_sigma = -deriv2;
        arma::mat B0_N_sigma = arma::inv_sympd(B0_N_inv_sigma);

        B0s_N_mu.col(k).subvec(0, dh0-1) = B0_N_mu;
        B0s_N_sigma.slice(k).submat(0, 0, dh0-1, dh0-1) = B0_N_sigma;

        B0v.col(k) = B0_.t() * (B0_mu-h0);

        // V_N_mu.col(k) = v;
        // V_N_sigma.slice(k) = -1/deriv2;
        //
        // double l_cmult = l_multinomial_const(X.col(k));
        //
        // int dh0_red = 1;
        // unsigned int index[] = {0,0};
        // int position = 0, kappa = 0;
        // arma::vec l_dens(order);
        //
        // double w;
        // do{
        //   v = v_n_mu + sqrt(v_n_sigma(0,0)) * uni_hermite(index[0],0);
        //   w = uni_hermite(index[0],dh0_red);
        //
        //
        //   h0_ = B0_mu + V * v;
        //   if(d1s(k)>0){
        //     h = B0s * h0_ + B1s * h1;
        //   }else{
        //     h = B0s * h0_;
        //   }
        //   h = h - h.max();
        //   arma::vec p = exp(h);
        //   l_dens(kappa) = w + l_dnormal_vec(h0_, B0_mu, inv_B0_sigma) -
        //     l_dnormal_vec(v, v_n_mu, inv_v_n_sigma) +
        //     arma::dot(log(p/accu(p)), X.col(k)) + l_cmult;
        //
        //   // Calculate next coordinate
        //   index[position]++;
        //   while(index[position] == order){
        //     index[position] = 0;
        //     position++;
        //     index[position]++;
        //   }
        //   position = 0;
        //   kappa++;
        // } while (index[(unsigned)dh0_red] == 0);
        // arma::vec dens = exp(l_dens-l_dens.max());
        //
        // double M0 = 0;
        // arma::vec M1 = arma::zeros(dh0_red);
        // arma::mat M2 = arma::zeros(dh0_red,dh0_red);
        // for(unsigned int i = 0; i <= dh0_red; i++) index[i] = 0;
        // position = 0;
        // kappa = 0;
        // do{
        //   v = v_n_mu + sqrt(v_n_sigma(0,0)) * uni_hermite(index[0],0);
        //
        //   M0 += dens(kappa);
        //   M1 += v * dens(kappa);
        //   M2 += v.t() * v * dens(kappa);
        //
        //   // Calculate next coordinate
        //   index[position]++;
        //   while(index[position] == order){
        //     index[position] = 0;
        //     position++;
        //     index[position]++;
        //   }
        //   position = 0;
        //   kappa++;
        // } while (index[(unsigned)dh0_red] == 0);
        // EV1.col(k) = M1/M0;
        // EV2.slice(k) = M2/M0;
        //
        // // cbind(B0,B1) %*% rbind(
        // //     cbind(V %*% cbind(I2$EV2[,,1]) %*% t(V), V %*% cbind(I2$EV1[,1]) %*% t(h1)),
        // //     cbind(h1 %*% cbind(I2$EV1[,1]) %*% t(V), h1 %*% t(h1))) %*% t(cbind(B0,B1))
        //
        //
        // arma::vec Eh0 = B0_mu + V * (M1/M0);
        // arma::mat Eh0h0 = B0_mu * B0_mu.t() + B0_mu * (V * M1/M0).t() + V * (M1/M0) * B0_mu.t() + V * (M2/M0) * V.t();
        //
        // if(d1s(k)>0){
        //   arma::mat B = arma::join_rows(B0s,B1s);
        //   clr_E1.col(k) = B * arma::join_cols(Eh0, h1);
        //   clr_E2.slice(k) = B * arma::join_cols(
        //     arma::join_rows(Eh0h0, Eh0 * h1.t()),
        //     arma::join_rows(h1 * Eh0.t(), h1 * h1.t())) * B.t();
        // }else{
        //   clr_E1.col(k) = B0s * Eh0;
        //   clr_E2.slice(k) = B0s * Eh0h0 * B0s.t();
        // }
      }
    }
    clr_mu_new = mean(clr_E1, 1);
    em_current_eps = norm(clr_mu-clr_mu_new, 1);

    clr_mu = clr_mu_new;
    arma::mat M = mean(clr_E2, 2);
    clr_sigma = M - clr_mu_new * clr_mu_new.t();
  }while(em_iter < em_max_iter & em_current_eps > em_eps);

  lrnm_res["d0"] = d0;
  lrnm_res["d1"] = d1;
  lrnm_res["d0s"] = d0s;
  lrnm_res["d1s"] = d1s;
  lrnm_res["I0"] = I0;
  lrnm_res["I1"] = I1;
  lrnm_res["B0"] = B0;
  lrnm_res["B1"] = B1;
  lrnm_res["P"] = P;
  lrnm_res["X"] = X;
  lrnm_res["clr_H"] = clr_E1;
  lrnm_res["clr_mu"] = clr_mu;
  lrnm_res["clr_sigma"] = clr_sigma;
  lrnm_res["B0v"] = B0v;

  lrnm_res["B0s_N_mu"] = B0s_N_mu;
  lrnm_res["B0s_N_sigma"] = B0s_N_sigma;
  lrnm_res["clr_E1"] = clr_E1;
  lrnm_res["clr_E2"] = clr_E2;
  lrnm_res["em_iter"] = em_iter;
  return(lrnm_res);
}


// [[Rcpp::export]]
List c_cond_lrnm_V_fit_hermite(arma::mat& X, arma::mat& C,
                               arma::mat& Vclr,
                               int order,
                               double em_eps, int em_max_iter){
  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;

  arma::vec d0 = arma::vec(n);
  arma::umat I0 = arma::umat(D, n);
  for(int k=0; k<n; k++){
    d0(k) = arma::accu(C.col(k) == 0);
    I0.col(k).head(d0(k)) = arma::find(C.col(k) == 0);
  }

  arma::vec d1 = d-d0;
  arma::vec D1 = 1+d1;
  arma::umat I1 = arma::umat(D, n);
  for(int k=0; k<n; k++){
    I1.col(k).head(D1(k)) = arma::find(C.col(k) > 0);
  }

  List lrnm_res = c_dm_fit(X, 0.001, 1000);
  arma::vec alpha = lrnm_res["alpha"];
  arma::mat P(X);
  for(unsigned int k=0;k<n;k++){
    arma::uvec uk(1);
    uk[0] = k;
    arma::vec a_p = X.col(k) + alpha;
    arma::vec a_p_1 = a_p / arma::accu(a_p);
    double res = arma::accu(a_p_1(I1.col(k).head(D1(k))));
    arma::vec a_p_0 = a_p(I0.col(k).head(d0(k)));
    a_p_0 = a_p_0 / arma::accu(a_p_0) * (1-res);
    P.col(k) = P.col(k) / arma::accu(P.col(k)) * res;
    P(I0.col(k).head(d0(k)), uk) = a_p_0;
  }


  arma::mat clr_H = log(P);
  arma::rowvec clr_m = arma::mean(clr_H, 0);
  for(int k=0; k<n; k++) clr_H.col(k) -= clr_m(k);

  arma::vec clr_mu = arma::mean(clr_H, 1);
  arma::mat clr_sigma = arma::cov(clr_H.t());

  arma::cube B1 = arma::zeros(D,d,n);
  arma::vec d1s = arma::zeros(n);
  arma::vec d1_unique = unique(d1);
  arma::cube ILR1 = arma::zeros(D,d,d1_unique.n_elem);

  arma::vec i1 = arma::zeros(D+1);
  for(int i=0;i<d1_unique.n_elem;i++){
    i1(d1_unique(i)) = i;
    if(d1_unique(i) > 0) ILR1.slice(i).submat(0,0,d1_unique(i),d1_unique(i)-1) = ilr_basis_default(d1_unique(i)+1);
  }
  for(int k=0;k<n;k++){
    if(d1(k) > 0){
      arma::uvec rows = I1.col(k).head(D1(k));
      arma::mat ilr_m = ILR1.slice(i1(d1(k))).head_rows(D1(k));
      B1.slice(k).rows(rows) = ilr_m;

      arma::mat B1s = B1.slice(k).head_cols(d1(k));
      arma::mat B1_sigma = B1s.t() * clr_sigma * B1s;

      arma::vec B1_eigval;
      arma::mat B1_eigvec;
      eig_sym(B1_eigval, B1_eigvec, B1_sigma);
      arma::mat B1_ = B1_eigvec.cols(arma::find(B1_eigval > 1e-3));
      d1s(k) = B1_.n_cols;
      B1.slice(k).head_cols(d1s(k)) = B1.slice(k).head_cols(d1(k)) * B1_;
    }
  }


  arma::cube B0 = arma::zeros(D,d,n);
  arma::vec d0s = arma::zeros(n);
  arma::vec d0_unique = unique(d0);
  arma::cube ILR0 = arma::zeros(D,d,d0_unique.n_elem);

  arma::vec i0 = arma::zeros(D+1);
  for(int i=0;i<d0_unique.n_elem;i++){
    i0(d0_unique(i)) = i;
    if(d0_unique(i)>1){
      ILR0.slice(i).submat(0,0,d0_unique(i)-1,d0_unique(i)-2) = ilr_basis_default(d0_unique(i));
    }
  }
  for(int k=0;k<n;k++){

    if(d0(k) > 0){
      if(d0(k) > 1){
        arma::uvec rows = I0.col(k).head(d0(k));
        arma::mat ilr_m = ILR0.slice(i0(d0(k))).head_rows(d0(k));
        B0.slice(k).rows(rows) = ilr_m;
      }
      double v0 = +1/d0(k) * std::sqrt(d0(k) * (D1(k))/D);
      for(int i=0; i<d0(k);i++){
        B0.slice(k)(I0(i,k),d0(k)-1) = v0;
      }
      double v1 = -1/(D1(k)) * std::sqrt(d0(k) * (D1(k))/D);
      for(int i=0; i<D1(k);i++){
        B0.slice(k)(I1(i,k),d0(k)-1) = v1;
      }
      arma::mat B0s = B0.slice(k).head_cols(d0(k));
      arma::mat B0_sigma = B0s.t() * clr_sigma * B0s;
      // Basis reduction
      if(d1(k)>0){
        arma::mat B1s = B1.slice(k).head_cols(d1s(k));
        arma::vec h1 = B1s.t() * clr_H.col(k);
        arma::mat clr_sigma_B1s = clr_sigma * B1s;
        arma::mat B1_0 = B0s.t() * clr_sigma_B1s * arma::inv_sympd(B1s.t() * clr_sigma_B1s);
        // B1_0.print("B1_0");
        // B0_sigma.print("B0_sigma");
        B0_sigma = B0_sigma - B1_0 * clr_sigma_B1s.t() * B0s;
        // B0_sigma.print("B0_sigma");
      }
      arma::vec B0_eigval;
      arma::mat B0_eigvec;
      eig_sym(B0_eigval, B0_eigvec, B0_sigma);
      arma::mat B0_ = B0_eigvec.cols(arma::find(B0_eigval > 1e-3));
      d0s(k) = B0_.n_cols;
      B0.slice(k).head_cols(d0s(k)) = B0.slice(k).head_cols(d0(k)) * B0_;
    }
  }

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));

  arma::mat B0v = arma::zeros(d,n);
  arma::mat V_N_mu(1,n);
  arma::cube V_N_sigma(1,1,n);

  arma::mat EV1(1,n);
  arma::cube EV2(1,1,n);

  arma::mat clr_E1(clr_H);
  arma::cube clr_E2(D,D,n);
  for(int k=0;k<n;k++) clr_E2.slice(k) = clr_E1.col(k) * clr_E1.col(k).t();

  ///
  int em_iter = 0;
  double em_current_eps;
  arma::vec clr_mu_new;
  do{
    em_iter++;
    for(int k=0;k<n;k++){
      if(d0s(k)>0){
        arma::mat B0s = B0.slice(k).head_cols(d0s(k));
        arma::mat B1s = B1.slice(k).head_cols(d1s(k));
        arma::vec clr_h = clr_E1.col(k);
        arma::vec h1 = B1s.t() * clr_h;
        arma::vec h0 = B0s.t() * clr_h;

        // (B1s.t() * clr_sigma * B1s).print("B1s.t() * clr_sigma * B1s");
        arma::mat inv_B1_sigma = arma::inv_sympd(B1s.t() * clr_sigma * B1s);
        // inv_B1_sigma.print("inv_B1");
        arma::mat B1_0 = B1s.t() * B0s;

        arma::vec B0_mu = B0s.t() * clr_mu + B1_0.t() * inv_B1_sigma * (h1 - B1s.t() * clr_mu);
        arma::mat B0_sigma = B0s.t() * clr_sigma * B0s - B1_0.t() * inv_B1_sigma * B1_0;
        arma::mat inv_B0_sigma = arma::inv_sympd(B0_sigma);


        arma::vec V = B0s.t() * Vclr.col(k);

        B0v.col(k).head(d0s(k)) = V;

        int current_iter = 0;
        int dh0 = 1;
        arma::vec deriv1(1), deriv2(1), step = arma::ones(1);
        arma::vec eps = {1e-5};
        int max_iter = 100;

        arma::vec v = V.t() * h0; // {0}; //
        arma::vec h0_(d0s(k));
        arma::vec h(D);
        do{

          current_iter++;
          h0_ = h0 + V * v;
          if(d1s(k)>0){
            h = B0s * h0_ + B1s * h1;
          }else{
            h = B0s * h0_;
          }
          arma::vec eBh = exp(h);
          double  w = arma::accu(eBh);

          arma::vec  wBi = B0s.t() * eBh;
          deriv1 = V.t() * (-inv_B0_sigma  * (h0_-B0_mu) + B0s.t() * X.col(k) - sum(X.col(k)) * wBi / w);

          arma::mat wBij(d0s(k), d0s(k));
          for(int i=0;i<d0s(k);i++){
            for(int j=0;j<d0s(k);j++){
              wBij(i,j) = arma::accu(B0s.col(i) % B0s.col(j) % eBh);
            }
          }

          deriv2 = V.t() * (-inv_B0_sigma - sum(X.col(k)) * ( -(wBi * wBi.t())/ (w*w) + wBij / w)) * V;
          step = deriv1/deriv2;
          v = v - 0.9 * step;

        }while( std::abs(step(0)) > eps(0) && current_iter < max_iter);
        arma::vec v_n_mu = v;
        arma::mat v_n_sigma = -1/deriv2;
        arma::mat inv_v_n_sigma = -deriv2;

        V_N_mu.col(k) = v;
        V_N_sigma.slice(k) = -1/deriv2;

        double l_cmult = l_multinomial_const(X.col(k));

        int dh0_red = 1;
        unsigned int index[] = {0,0};
        int position = 0, kappa = 0;
        arma::vec l_dens(order);

        double w;
        do{
          v = v_n_mu + sqrt(v_n_sigma(0,0)) * uni_hermite(index[0],0);
          w = uni_hermite(index[0],dh0_red);


          h0_ = h0 + V * v;
          if(d1s(k)>0){
            h = B0s * h0_ + B1s * h1;
          }else{
            h = B0s * h0_;
          }
          h = h - h.max();
          arma::vec p = exp(h);
          l_dens(kappa) = w + l_dnormal_vec(h0_, B0_mu, inv_B0_sigma) -
            l_dnormal_vec(v, v_n_mu, inv_v_n_sigma) +
            arma::dot(log(p/accu(p)), X.col(k)) + l_cmult;

          // Calculate next coordinate
          index[position]++;
          while(index[position] == order){
            index[position] = 0;
            position++;
            index[position]++;
          }
          position = 0;
          kappa++;
        } while (index[(unsigned)dh0_red] == 0);
        arma::vec dens = exp(l_dens-l_dens.max());

        double M0 = 0;
        arma::vec M1 = arma::zeros(dh0_red);
        arma::mat M2 = arma::zeros(dh0_red,dh0_red);
        for(unsigned int i = 0; i <= dh0_red; i++) index[i] = 0;
        position = 0;
        kappa = 0;
        do{
          v = v_n_mu + sqrt(v_n_sigma(0,0)) * uni_hermite(index[0],0);

          M0 += dens(kappa);
          M1 += v * dens(kappa);
          M2 += v.t() * v * dens(kappa);

          // Calculate next coordinate
          index[position]++;
          while(index[position] == order){
            index[position] = 0;
            position++;
            index[position]++;
          }
          position = 0;
          kappa++;
        } while (index[(unsigned)dh0_red] == 0);
        EV1.col(k) = M1/M0;
        EV2.slice(k) = M2/M0;

        // cbind(B0,B1) %*% rbind(
        //     cbind(V %*% cbind(I2$EV2[,,1]) %*% t(V), V %*% cbind(I2$EV1[,1]) %*% t(h1)),
        //     cbind(h1 %*% cbind(I2$EV1[,1]) %*% t(V), h1 %*% t(h1))) %*% t(cbind(B0,B1))


        arma::vec Eh0 = h0 + V * (M1/M0);
        arma::mat Eh0h0 = h0 * h0.t() + h0 * (V * M1/M0).t() + V * (M1/M0) * h0.t() + V * (M2/M0) * V.t();

        if(d1s(k)>0){
          arma::mat B = arma::join_rows(B0s,B1s);
          clr_E1.col(k) = B * arma::join_cols(Eh0, h1);
          clr_E2.slice(k) = B * arma::join_cols(
            arma::join_rows(Eh0h0, Eh0 * h1.t()),
            arma::join_rows(h1 * Eh0.t(), h1 * h1.t())) * B.t();
        }else{
          clr_E1.col(k) = B0s * Eh0;
          clr_E2.slice(k) = B0s * Eh0h0 * B0s.t();
        }
      }
    }
    clr_mu_new = mean(clr_E1, 1);
    em_current_eps = norm(clr_mu-clr_mu_new, 1);

    clr_mu = clr_mu_new;
    arma::mat M = mean(clr_E2, 2);
    clr_sigma = M - clr_mu_new * clr_mu_new.t();
  }while(em_iter < em_max_iter & em_current_eps > em_eps);

  lrnm_res["d0"] = d0;
  lrnm_res["d1"] = d1;
  lrnm_res["d0s"] = d0s;
  lrnm_res["d1s"] = d1s;
  lrnm_res["I0"] = I0;
  lrnm_res["I1"] = I1;
  lrnm_res["B0"] = B0;
  lrnm_res["B1"] = B1;
  lrnm_res["P"] = P;
  lrnm_res["X"] = X;
  lrnm_res["clr_H"] = clr_E1;
  lrnm_res["clr_mu"] = clr_mu;
  lrnm_res["clr_sigma"] = clr_sigma;
  lrnm_res["B0v"] = B0v;

  lrnm_res["V_N_mu"] = V_N_mu;
  lrnm_res["V_N_sigma"] = V_N_sigma;
  lrnm_res["clr_E1"] = clr_E1;
  lrnm_res["clr_E2"] = clr_E2;
  lrnm_res["em_iter"] = em_iter;
  return(lrnm_res);
}

