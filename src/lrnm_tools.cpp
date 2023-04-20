// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include "dm.h"
#include "coda_base.h"
#include "lrnm_utils.h"
#include "hermite.h"

using namespace Rcpp;

// [[Rcpp::export]]
List c_cond_lrnm_init(arma::mat& X, arma::mat& C, int d1max = 0, int d0max = 0){
  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;

  if(d0max == 0) d0max = d;
  if(d1max == 0) d1max = d;

  d0max = std::min(d0max, n-1);
  d1max = std::min(d1max, n-1);

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

  List lrnm_ini = c_dm_fit(X, 0.001, 1000);
  arma::vec alpha = lrnm_ini["alpha"];
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
    ILR1.slice(i).submat(0,0,d1_unique(i),d1_unique(i)-1) = ilr_basis_default(d1_unique(i)+1);
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
      if(B1_.n_cols > d1max){
        B1_ = B1_.tail_cols(d1max);
      }
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

      if(B0_.n_cols > d0max){
        B0_ = B0_.tail_cols(d0max);
      }
      d0s(k) = B0_.n_cols;
      B0.slice(k).head_cols(d0s(k)) = B0.slice(k).head_cols(d0(k)) * B0_;
    }
  }

  lrnm_ini["d0"] = d0;
  lrnm_ini["d1"] = d1;
  lrnm_ini["d0s"] = d0s;
  lrnm_ini["d1s"] = d1s;
  lrnm_ini["I0"] = I0;
  lrnm_ini["I1"] = I1;
  lrnm_ini["B0"] = B0;
  lrnm_ini["B1"] = B1;
  lrnm_ini["P"] = P;
  lrnm_ini["X"] = X;
  lrnm_ini["clr_H"] = clr_H;
  lrnm_ini["clr_mu"] = clr_mu;
  lrnm_ini["clr_sigma"] = clr_sigma;
  return(lrnm_ini);
}

// [[Rcpp::export]]
List c_cond_lrnm_laplace(List lrnm_model){
  arma::mat X = lrnm_model["X"];
  arma::mat clr_H = lrnm_model["clr_H"];
  arma::vec clr_mu = lrnm_model["clr_mu"];
  arma::mat clr_sigma = lrnm_model["clr_sigma"];
  arma::cube B0 = lrnm_model["B0"];
  arma::cube B1 = lrnm_model["B1"];
  arma::vec d0s = lrnm_model["d0s"];
  arma::vec d1s = lrnm_model["d1s"];

  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;

  arma::mat B0_N_mu = arma::zeros(d, n);
  arma::cube B0_N_sigma = arma::zeros(d, d, n);
  for(int k=0;k<n;k++){
    if(d0s(k)>0){
      arma::mat T0 = B0.slice(k).head_cols(d0s(k)).t();
      arma::mat T1 = B1.slice(k).head_cols(d1s(k)).t();

      arma::vec B0_mu = T0 * clr_mu;
      arma::mat B0_sigma = T0 * clr_sigma * T0.t();
      arma::mat inv_B0_sigma = arma::inv_sympd(B0_sigma);

      int current_iter = 0;
      int dh0 = d0s(k);


      arma::vec h0 = T0 * clr_H.col(k);
      arma::vec h1 = T1 * clr_H.col(k);

      arma::vec deriv1(dh0);
      arma::mat deriv2(dh0, dh0);
      arma::vec step = arma::ones(dh0);
      double eps = 1e-5;
      int max_iter = 100;

      // if(k==2) T0.print("T0");
      do{

        current_iter++;

        arma::vec h = T0.t() * h0;
        if(d1s(k)>0){
          h += T1.t() * h1;
        }
        arma::vec eBh = exp(h);

        double  w = arma::accu(eBh);

        arma::vec  wBi = T0 * eBh;
        deriv1 = -inv_B0_sigma  * (h0-B0_mu) + T0 * X.col(k) - sum(X.col(k)) * wBi / w;

        arma::mat wBij(dh0,dh0);
        for(int i=0; i<dh0; i++){
          for(int j=0; j<dh0; j++){
            wBij(i,j) = arma::accu(T0.row(i).t() % T0.row(j).t() % eBh);
          }
        }
        deriv2 = -inv_B0_sigma - sum(X.col(k)) * ( -(wBi * wBi.t())/ (w*w) + wBij / w);

        step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
        h0 = h0 - 0.9 * step;

      }while( norm(step, 2) > eps && current_iter < max_iter);
      arma::mat B0_N_inv_sigma = -deriv2;
      B0_N_mu.col(k).subvec(0, d0s(k)-1) = h0;
      B0_N_sigma.slice(k).submat(0, 0, d0s(k)-1, d0s(k)-1) = arma::inv_sympd(B0_N_inv_sigma);
    }
  }
  lrnm_model["B0_N_mu"] = B0_N_mu;
  lrnm_model["B0_N_sigma"] = B0_N_sigma;
  return(lrnm_model);

}

// [[Rcpp::export]]
List c_cond_lrnm_pc_expected_hermite(List lrnm_model, int d0red = 1, int order = 10){
  arma::mat X = lrnm_model["X"];
  arma::mat clr_H = lrnm_model["clr_H"];
  arma::vec clr_mu = lrnm_model["clr_mu"];
  arma::mat clr_sigma = lrnm_model["clr_sigma"];
  arma::cube B0 = lrnm_model["B0"];
  arma::cube B1 = lrnm_model["B1"];
  arma::vec d0s = lrnm_model["d0s"];
  arma::vec d1s = lrnm_model["d1s"];
  arma::mat B0_N_mu = lrnm_model["B0_N_mu"];
  arma::cube B0_N_sigma = lrnm_model["B0_N_sigma"];

  int n = X.n_cols;
  int D = X.n_rows;
  int d = D - 1;

  if(d0red == 0) d0red = d;

  arma::mat clr_E1(D,n);
  arma::cube clr_E2(D,D,n);

  arma::mat uni_hermite = hermite(order);
  uni_hermite.col(1) = log(uni_hermite.col(1));
  for(int k=0;k<n;k++){
    arma::vec B0_N_mu_k = B0_N_mu.col(k).subvec(0, d0s(k)-1);
    arma::mat B0_N_sigma_k = B0_N_sigma.slice(k).submat(0, 0, d0s(k)-1, d0s(k)-1);
    arma::vec B0s_eigval;
    arma::mat B0s_eigvec;
    eig_sym(B0s_eigval, B0s_eigvec, B0_N_sigma_k);

    int dh0_red = std::min((int)d0s(k), d0red);

    arma::mat T0sub = B0s_eigvec.tail_cols(dh0_red).t();
    arma::mat T0 = T0sub * B0.slice(k).head_cols(d0s(k)).t();
    arma::mat T1 = B1.slice(k).head_cols(d1s(k)).t();
    // B0s_eigval.t().print("B0s_eigval");
    // B0s_eigvec.print("B0s_eigvec");
    // arma::mat rotation = arma::diagmat(flipud(sqrt(B0s_eigval.tail(dh0_red))));
    arma::vec B0_red_mu = T0 * clr_mu;
    arma::mat B0_red_sigma = T0 * clr_sigma * T0.t();
    arma::mat B0_red_inv_sigma = arma::inv_sympd(B0_red_sigma);

    arma::vec B0_red_N_mu = T0sub * B0_N_mu_k;
    arma::mat B0_red_N_sigma = T0sub * B0_N_sigma_k * T0sub.t();
    arma::mat B0_red_N_inv_sigma = T0sub * arma::inv_sympd(B0_N_sigma_k) * T0sub.t();

    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, B0_red_N_sigma);
    // eigval.t().print("eigval");
    // eigvec.print("eigvec");
    arma::mat rotation = fliplr(eigvec) * arma::diagmat(flipud(sqrt(eigval)));

    double l_cmult = l_multinomial_const(X.col(k));

    arma::vec h0_red = arma::zeros(dh0_red);
    arma::vec h1 = T1 * clr_H.col(k);

    unsigned int index[(unsigned)dh0_red+1];
    for(unsigned int i = 0; i <= dh0_red; i++) index[i] = 0;
    int position = 0, kappa = 0;
    arma::vec l_dens(std::pow(order, dh0_red));
    double l_dens_max = 0;

    do{
      double w = 0;
      for(unsigned int i = 0; i < dh0_red; i++){
        h0_red(i) = uni_hermite(index[i],0);
        w += uni_hermite(index[i],1);
      }
      h0_red = B0_red_N_mu + rotation * h0_red;

      arma::vec h = T0.t() * h0_red;
      if(d1s(k)>0){
        h += T1.t() * h1;
      }

      arma::vec Bh = h;
      Bh = Bh - Bh.max();
      arma::vec p = exp(Bh);
      l_dens(kappa) = w + l_dnormal_vec(h0_red, B0_red_mu, B0_red_inv_sigma) -
        l_dnormal_vec(h0_red, B0_red_N_mu, B0_red_N_inv_sigma) +
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
      for(unsigned int i = 0; i < dh0_red; i++){
        h0_red(i) = uni_hermite(index[i],0);
      }
      h0_red = B0_red_N_mu + rotation * h0_red;

      M0 += dens(kappa);
      M1 += h0_red * dens(kappa);
      M2 += h0_red * h0_red.t() * dens(kappa);

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
    clr_E1.col(k) = T0.t() * (M1/M0);
    clr_E2.slice(k) = T0.t() * (M2/M0) * T0;
    arma::vec clr_h1 = T1.t() * h1;
    if(d1s(k)>0){
      clr_E1.col(k) += clr_h1;
      clr_E2.slice(k) += clr_h1 * clr_h1.t();
    }
  }
  lrnm_model["clr_E1"] = clr_E1;
  lrnm_model["clr_E2"] = clr_E2;
  return(lrnm_model);

}

// [[Rcpp::export]]
List c_cond_lrnm_fixed_point_expected_hermite(List lrnm_model, int order = 10){
 arma::mat X = lrnm_model["X"];
 arma::mat clr_H = lrnm_model["clr_H"];
 arma::vec clr_mu = lrnm_model["clr_mu"];
 arma::mat clr_sigma = lrnm_model["clr_sigma"];
 arma::cube B0 = lrnm_model["B0"];
 arma::cube B1 = lrnm_model["B1"];
 arma::vec d0s = lrnm_model["d0s"];
 arma::vec d1s = lrnm_model["d1s"];

 int n = X.n_cols;
 int D = X.n_rows;
 int d = D - 1;

 arma::mat V_N_mu(1,n);
 arma::cube V_N_sigma(1,1,n);

 arma::mat EV1(1,n);
 arma::cube EV2(1,1,n);

 arma::mat clr_E1(D,n);
 arma::cube clr_E2(D,D,n);
 arma::mat B0v = arma::zeros(d,n);


 arma::mat uni_hermite = hermite(order);
 uni_hermite.col(1) = log(uni_hermite.col(1));
 for(int k=0;k<n;k++){
   if(d0s(k)>0){
     arma::mat B0s = B0.slice(k).head_cols(d0s(k));
     arma::mat B1s = B1.slice(k).head_cols(d1s(k));
     arma::vec clr_h = clr_H.col(k);
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

     arma::vec v = {0};
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
   }else{ // d0s(k)==0
     arma::mat B1s = B1.slice(k).head_cols(d1s(k));
     arma::vec h1 = B1s.t() * clr_H.col(k);

     clr_E1.col(k) = B1s * h1;
     clr_E2.slice(k) = B1s * h1 * h1.t() * B1s.t();
   }
 }

 lrnm_model["B0v"] = B0v;

 lrnm_model["V_N_mu"] = V_N_mu;
 lrnm_model["V_N_sigma"] = V_N_sigma;

 lrnm_model["EV1"] = EV1;
 lrnm_model["EV2"] = EV2;

 lrnm_model["clr_E1"] = clr_E1;
 lrnm_model["clr_E2"] = clr_E2;
 return(lrnm_model);

}

