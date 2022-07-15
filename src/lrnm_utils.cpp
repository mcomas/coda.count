// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

using namespace Rcpp;

const double log2pi = std::log(M_2PI);

// [[Rcpp::export]]
double l_dnormal_vec(arma::vec h, arma::vec mu, arma::mat inv_sigma){
  int k = h.n_elem;

  double log_det_val;
  double sign;
  log_det(log_det_val, sign, inv_sigma);
  double norm_const = -0.5 * k * log2pi + 0.5 * log_det_val;
  double log_norm = 0;

  arma::vec x;
  x = h - mu;
  x = x.t() * inv_sigma * x;
  log_norm = -0.5 * x(0);

  return(norm_const + log_norm);
}

// [[Rcpp::export]]
double l_dnormal_prop_vec(arma::vec h, arma::vec mu, arma::mat inv_sigma){
  double log_norm = 0;
  arma::vec x;
  x = h - mu;
  x = x.t() * inv_sigma * x;
  log_norm = -0.5 * x(0);

  return(log_norm);
}


double l_multinomial_const(arma::vec x){
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

double l_multinomial(arma::vec x, arma::vec p, double lconst){
  //double lconst = lpmultinomial_const(x);
  return( lconst + arma::dot(log(p),x) );
}

// [[Rcpp::export]]
double l_lrnm_join_no_constant_vec(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){

  arma::vec Bh = Binv * h;
  double lmult = arma::accu( x % (Bh - log(arma::accu(exp(Bh))) ) );
  double lnormal = l_dnormal_vec(h, mu, inv_sigma);

  return(lmult + lnormal);
}

// double neg_l_lrnm_join_no_constant_vec(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){
//
//   arma::vec Bh = Binv * h;
//   double lmult = arma::accu( x % (Bh - log(arma::accu(exp(Bh))) ) );
//   double lnormal = l_dnormal_vec(h, mu, inv_sigma);
//
//   return(-lmult - lnormal);
// }

// [[Rcpp::export]]
double l_lrnm_join_vec(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){

  arma::vec Bh = Binv * h;
  double lconst = l_multinomial_const(x);
  double lmult = arma::accu( x % (Bh - log(arma::accu(exp(Bh))) ) );
  double lnormal = l_dnormal_vec(h, mu, inv_sigma);

  return(lconst + lmult + lnormal);
}

//' @export
// [[Rcpp::export]]
arma::vec l_lrnm_join_d1(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){
  int k = h.size();
  arma::vec deriv(k);

  arma::vec eBh = exp(Binv * h);

  double  w = arma::accu(eBh);
  arma::vec  wBi = Binv.t() * eBh;
  return(-inv_sigma * (h-mu) + Binv.t() * x - sum(x) * wBi / w);
}

//' @export
// [[Rcpp::export]]
arma::vec l_lrnm_cond_join_d1(arma::vec h1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv){
  int k = h1.size();
  arma::vec deriv(k);


  arma::vec h = join_cols(h1,h2);
  arma::vec eBh = exp(Binv * h);

  double  w = arma::accu(eBh);

  arma::vec  wBi = Binv.head_cols(k).t() * eBh;

  return(-inv_sigma * (h1-mu) + Binv.head_cols(k).t() * x - sum(x) * wBi / w);
}

//' @export
// [[Rcpp::export]]
arma::vec l_indep_lrnm_cond_join_d1(arma::vec h1, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv){
  arma::vec h = join_cols(h1,h2);
  arma::vec eBh = exp(Binv * h);

  double  w = arma::accu(eBh);

  arma::vec  wBi = Binv.col(ih1).t() * eBh;

  return(-inv_sigma * (h1(ih1)-mu) + Binv.col(ih1).t() * x - sum(x) * wBi / w);
}

//' @export
// [[Rcpp::export]]
arma::vec l_lrnm_cond1_join_d1(arma::vec h, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){

  arma::vec eBh = exp(Binv * h);

  double  w = arma::accu(eBh);

  arma::vec  wBi = Binv.col(ih1).t() * eBh;

  return(-inv_sigma * (h(ih1)-mu) + Binv.col(ih1).t() * x - sum(x) * wBi / w);
}

// arma::vec neg_l_lrnm_join_d1(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){
//   int k = h.size();
//   arma::vec deriv(k);
//
//   arma::vec eBh = exp(Binv * h);
//
//   double  w = arma::accu(eBh);
//   arma::vec  wBi = Binv.t() * eBh;
//   return(inv_sigma * (h-mu) - Binv.t() * x + sum(x) * wBi / w);
// }

//' @export
// [[Rcpp::export]]
arma::mat l_lrnm_join_d2(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){
  int k = h.size();
  arma::mat deriv(k,k);

  arma::vec eBh = exp(Binv * h);
  double  w = arma::accu(eBh);
  arma::vec wBi = Binv.t() * eBh;
  arma::mat wBij(k,k);
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      wBij(i,j) = arma::accu(Binv.col(i) % Binv.col(j) % eBh);
    }
  }
  return(-inv_sigma - sum(x) * ( -(wBi * wBi.t())/ (w*w) + wBij / w));
}

//' @export
// [[Rcpp::export]]
arma::mat l_lrnm_cond_join_d2(arma::vec h1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv){
  int k = h1.size();
  arma::mat deriv(k,k);

  arma::vec h = join_cols(h1,h2);
  arma::vec eBh = exp(Binv * h);
  double  w = arma::accu(eBh);

  arma::vec wBi = Binv.head_cols(k).t() * eBh;
  arma::mat wBij(k,k);
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      wBij(i,j) = arma::accu(Binv.col(i) % Binv.col(j) % eBh);
    }
  }
  return(-inv_sigma - sum(x) * ( -(wBi * wBi.t())/ (w*w) + wBij / w));
}

//' @export
// [[Rcpp::export]]
arma::mat l_indep_lrnm_cond_join_d2(arma::vec h1, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv){


  arma::vec h = join_cols(h1,h2);
  arma::vec eBh = exp(Binv * h);
  double  w = arma::accu(eBh);

  arma::vec wBi = Binv.col(ih1).t() * eBh;

  double wBij = arma::accu(Binv.col(ih1) % Binv.col(ih1) % eBh);

  return(-inv_sigma - sum(x) * ( -(wBi * wBi.t())/ (w*w) + wBij / w));
}

//' @export
// [[Rcpp::export]]
arma::mat l_lrnm_cond1_join_d2(arma::vec h, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv){

  arma::vec eBh = exp(Binv * h);
  double  w = arma::accu(eBh);

  arma::vec wBi = Binv.col(ih1).t() * eBh;

  double wBij = arma::accu(Binv.col(ih1) % Binv.col(ih1) % eBh);

  return(-inv_sigma - sum(x) * ( -(wBi * wBi.t())/ (w*w) + wBij / w));
}

//' @export
// [[Rcpp::export]]
arma::vec l_lrnm_join_maximum(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv,
                              double eps = 1e-8, int max_iter = 1000){

  int k = Binv.n_cols;
  arma::vec h;
  if(x.min() > 0 & x.max() > 5){
    h = arma::pinv(Binv) * log(x);
  }else if(accu(x) > 100){
    h = arma::pinv(Binv) * log(x+1);
  }else{
    h = arma::vec(mu);
    h = 0.5 * arma::pinv(Binv) * log(x+1) + 0.5 * arma::vec(mu);
  }

  arma::vec deriv1(k);
  arma::mat deriv2(k,k);
  arma::vec step = arma::zeros<arma::vec>(k);

  int current_iter = 0;
  do{

    current_iter++;

    deriv1 = l_lrnm_join_d1(h, x, mu, inv_sigma, Binv);
    deriv2 = l_lrnm_join_d2(h, x, mu, inv_sigma, Binv);

    // Rcpp::Rcout << h << std::endl << deriv1 << std::endl << deriv2 << std::endl;
    step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
    h = h - 0.9 * step;
  }while( norm(step, 2) > eps && current_iter < max_iter);

  return h;
}

//' @export
// [[Rcpp::export]]
arma::vec l_lrnm_cond_join_maximum(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv,
                                   double eps = 1e-8, int max_iter = 1000){

  int k = mu.size();
  arma::vec h;
  if(x.min() > 0 & x.max() > 5){
    h = arma::pinv(Binv.head_cols(k)) * log(x);
  }else if(accu(x) > 100){
    h = arma::pinv(Binv.head_cols(k)) * log(x+1);
  }else{
    h = arma::vec(mu);
    h = 0.5 * arma::pinv(Binv.head_cols(k)) * log(x+1) + 0.5 * arma::vec(mu);
  }

  arma::vec deriv1(k);
  arma::mat deriv2(k,k);
  arma::vec step = arma::zeros<arma::vec>(k);

  int current_iter = 0;
  do{

    current_iter++;

    deriv1 = l_lrnm_cond_join_d1(h, x, mu, inv_sigma, h2, Binv);
    deriv2 = l_lrnm_cond_join_d2(h, x, mu, inv_sigma, h2, Binv);

    // Rcpp::Rcout << h << std::endl << deriv1 << std::endl << deriv2 << std::endl;
    step = arma::solve(deriv2, deriv1, arma::solve_opts::fast);
    h = h - 0.9 * step;
  }while( norm(step, 2) > eps && current_iter < max_iter);

  return h;
}

//' @export
// [[Rcpp::export]]
arma::vec l_indep_lrnm_cond_join_maximum(arma::vec h1, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv,
                                         double eps = 1e-8, int max_iter = 1000){



  arma::vec deriv1(1);
  arma::mat deriv2(1,1);
  double step = 0;
  arma::vec h(h1);
  int current_iter = 0;
  do{

    current_iter++;
    // Rcpp::Rcout << h << std::endl;
    deriv1 = l_indep_lrnm_cond_join_d1(h, ih1, x, mu, inv_sigma, h2, Binv);
    deriv2 = l_indep_lrnm_cond_join_d2(h, ih1, x, mu, inv_sigma, h2, Binv);

    // Rcpp::Rcout << h << std::endl << deriv1 << std::endl << deriv2 << std::endl;
    step = deriv1(0) / deriv2(0,0);
    h(ih1) = h(ih1) - 0.9 * step;
  }while( step > eps && current_iter < max_iter);

  return h;
}

//' @export
// [[Rcpp::export]]
double l_lrnm_cond1_join_maximum(arma::vec h, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv,
                                    double eps = 1e-8, int max_iter = 1000){



  arma::vec deriv1(1);
  arma::mat deriv2(1,1);
  double step = 0;
  int current_iter = 0;
  do{

    current_iter++;
    // Rcpp::Rcout << h << std::endl;
    deriv1 = l_lrnm_cond1_join_d1(h, ih1, x, mu, inv_sigma, Binv);
    deriv2 = l_lrnm_cond1_join_d2(h, ih1, x, mu, inv_sigma, Binv);

    // Rcpp::Rcout << h << std::endl << deriv1 << std::endl << deriv2 << std::endl;
    step = deriv1(0) / deriv2(0,0);
    h(ih1) = h(ih1) - 0.9 * step;
  }while( step > eps && current_iter < max_iter);

  return h(ih1);
}

// //' @export
// // [[Rcpp::export]]
// arma::vec l_lrnm_cond1_join_maximum(unsigned int ih1, arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv,
//                                     double eps = 1e-8, int max_iter = 1000){
//
//
//
//   arma::vec deriv1(1);
//   arma::mat deriv2(1,1);
//   double step = 0;
//   arma::vec h(h1);
//   int current_iter = 0;
//   do{
//
//     current_iter++;
//     // Rcpp::Rcout << h << std::endl;
//     deriv1 = l_indep_lrnm_cond_join_d1(h, ih1, x, mu, inv_sigma, h2, Binv);
//     deriv2 = l_indep_lrnm_cond_join_d2(h, ih1, x, mu, inv_sigma, h2, Binv);
//
//     // Rcpp::Rcout << h << std::endl << deriv1 << std::endl << deriv2 << std::endl;
//     step = deriv1(0) / deriv2(0,0);
//     h(ih1) = h(ih1) - 0.9 * step;
//   }while( step > eps && current_iter < max_iter);
//
//   return h;
// }
