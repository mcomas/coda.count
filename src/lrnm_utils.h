#include <RcppArmadillo.h>

#ifndef LRNM_UTILS_H
#define LRNM_UTILS_H



double l_lrnm_join_vec(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv);
double l_lrnm_join_no_constant_vec(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv);
arma::vec l_lrnm_join_d1(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv);
arma::mat l_lrnm_join_d2(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv);
arma::vec l_lrnm_join_maximum(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv,
                              double eps, int max_iter);


double l_dnormal_vec(arma::vec h, arma::vec mu, arma::mat& inv_sigma);
double l_multinomial_const(arma::vec x);
double l_multinomial(arma::vec x, arma::vec p, double lconst);

class LRNM_join{
public:
  // Construct the object with the given the design matrix and responses.
  LRNM_join(arma::vec& x, arma::vec& mu, arma::mat& inv_sigma, arma::mat& Binv) :
  x(x), mu(mu), inv_sigma(inv_sigma), Binv(Binv) { }
  void set_x(arma::vec& x_){
    x = x_;
  }
  void set_mu(arma::vec& mu_){
    mu = mu_;
  }
  // Return the objective function for model parameters beta.
  double Evaluate(const arma::mat& h)
  {
    return(-l_lrnm_join_no_constant_vec(h, x, mu, inv_sigma, Binv));
  }

  // Compute the gradient for model parameters beta
  void Gradient(const arma::mat& h, arma::mat& g)
  {
    g = -l_lrnm_join_d1(h, x, mu, inv_sigma, Binv);
  }

private:
  arma::vec& x;
  arma::vec& mu;
  arma::mat& inv_sigma;
  arma::mat& Binv;
};

arma::vec l_lrnm_cond_join_d1(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv);
arma::mat l_lrnm_cond_join_d2(arma::vec h, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv);
arma::vec l_lrnm_cond_join_maximum(arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv,
                                   double eps, int max_iter);

// arma::vec l_indep_lrnm_cond_join_d1(arma::vec h1, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv);
// arma::vec l_lrnm_cond1_join_d1(arma::vec h, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv);
//

// arma::vec l_indep_lrnm_cond_join_d2(arma::vec h1, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv);
// arma::vec l_lrnm_cond1_join_d2(arma::vec h, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv);
//

// arma::vec l_indep_lrnm_cond_join_maximum(arma::vec h1, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::vec h2, arma::mat &Binv,
//                                          double eps = 0.000001, int max_iter = 100);
// arma::vec l_lrnm_cond1_join_maximum(arma::vec h1, unsigned int ih1, arma::vec x, arma::vec mu, arma::mat &inv_sigma, arma::mat &Binv,
//                                          double eps = 0.000001, int max_iter = 100);

// double l_dnormal_prop_vec(arma::vec h, arma::vec mu, arma::mat inv_sigma);
#endif
