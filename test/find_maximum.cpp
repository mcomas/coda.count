#include <Rcpp.h>
#include "Eigen/Dense"
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include "cppoptlib/solver/bfgssolver.h"
#include "cppoptlib/solver/newtondescentsolver.h"

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

template <typename MatrixType>
inline typename MatrixType::Scalar logdet(const MatrixType& M, bool use_cholesky = false) {
  using namespace Eigen;
  using std::log;
  typedef typename MatrixType::Scalar Scalar;
  Scalar ld = 0;
  if (use_cholesky) {
    LLT<Eigen::Matrix<Scalar,Dynamic,Dynamic>> chol(M);
    auto& U = chol.matrixL();
    for (unsigned i = 0; i < M.rows(); ++i)
      ld += log(U(i,i));
    ld *= 2;
  } else {
    PartialPivLU<Eigen::Matrix<Scalar,Dynamic,Dynamic>> lu(M);
    auto& LU = lu.matrixLU();
    Scalar c = lu.permutationP().determinant(); // -1 or 1
    for (unsigned i = 0; i < LU.rows(); ++i) {
      const auto& lii = LU(i,i);
      if (lii < Scalar(0)) c *= -1;
      ld += log(abs(lii));
    }
    ld += log(c);
  }
  return ld;
}

template<typename T>
class ConditionalMode : public cppoptlib::Problem<T> {
public:
  using typename cppoptlib::Problem<T>::TVector;
  using typename cppoptlib::Problem<T>::THessian;
  using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

protected:
  const TVector x;
  const TVector mu;
  const MatrixType inv_sigma;
  const double log2pi = std::log(2.0 * M_PI);

public:
  ConditionalMode(const TVector &x_,
                  const TVector &mu_,
                  const MatrixType &sigma_) : x(x_), mu(mu_), inv_sigma(sigma_.llt().solve(Eigen::MatrixXd::Identity(sigma_.rows(), sigma_.rows()))) {}
  T value(const TVector &a) {
    int k = a.size();
    int K = k + 1;

    double norm_const = -0.5 * k * log2pi - 0.5 * logdet(inv_sigma, true);
    MatrixType log_norm =  -0.5 * (a-mu).transpose() * inv_sigma * (a-mu);

    double norm = norm_const + log_norm(0,0);

    double x_total = 0;
    for(int i = 0; i < K; i++) x_total += x[i];

    double x_total_fact = log(x_total);
    for(int l = 1; l < x_total; l++) x_total_fact += log(l);

    double x_parts_fact = 0;
    for(int i = 0; i <= K; i++){
      for(int l = 1; l <= x[i]; l++){
        x_parts_fact += log(l);
      }
    }

    double mult_const = 0;
    mult_const += (x_total_fact - x_parts_fact);

    double kappa = 1;
    for(int i = 0; i < k; i++) kappa += exp(a[i]);
    double multinom = -x[k] * log(kappa);
    for(int i = 0; i < k; i++) multinom += x[i] * ( a[i] - log(kappa));

    double mult = mult_const + multinom;
    return(-norm - mult);
  }
  void gradient(const TVector &a, TVector &grad){
    int k = a.size();
    double kappa = 1;
    for(int i = 0; i < k; i++) kappa += exp(a(i));

    for(int I=0; I < k; I++){
      Eigen::MatrixXd log_norm = -(a-mu).transpose() * inv_sigma.col(I);
      double mult = 0;
      mult += x(I) * (kappa-exp(a(I))) / kappa;
      for(int i = 0; i < I; i++) mult += x(i) * (-exp(a(I))) / kappa;
      for(int i = I+1; i < k+1; i++) mult += x(i) * (-exp(a(I))) / kappa;
      grad(I) = -log_norm(0) - mult;
    }
  }

  void hessian(const TVector &a, THessian &hessian) {
    int k = a.size();
    double kappa = 1;
    for(int i = 0; i < k; i++) kappa += exp(a(i));

    for(int I = 0; I < k; I++){
      for(int J = 0; J < k; J++){
        hessian(I,J) = inv_sigma(I, J);
        if(I == J){
          for(int i = 0; i < k + 1; i++) hessian(I,J) += x(i) * exp(a(I)) * (kappa-exp(a(I))) / (kappa*kappa);
        }else{
          for(int i = 0; i < k + 1; i++) hessian(I,J) -= x(i) * exp(a(I) + a(J)) / (kappa*kappa);
        }
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector example_BfgsSolver(NumericVector x, NumericVector mu, NumericMatrix sigma) {
  int k = mu.size();
  Eigen::VectorXd x_(x.size());
  Eigen::VectorXd mu_(mu.size());
  Eigen::MatrixXd sigma_(sigma.nrow(), sigma.ncol());
  for(int i=0; i<x.size(); i++) x_(i) = x(i);
  for(int i=0; i<mu.size(); i++) mu_(i) = mu(i);
  for(int i=0; i<sigma.nrow(); i++) for(int j=0; j<sigma.ncol(); j++) sigma_(i, j) = sigma(i, j);

  typedef ConditionalMode<double> TConditionalMode;
  typedef typename TConditionalMode::TVector TVector;
  typedef typename TConditionalMode::MatrixType MatrixType;

  TConditionalMode f(x_, mu_, sigma_);
  cppoptlib::BfgsSolver<TConditionalMode> solver;
  Eigen::VectorXd x0(mu.size());
  for(int i=0; i<mu.size(); i++){
    x0(i) = log(x(i)+0.5) - log(x(k)+0.5);
  }
  solver.minimize(f, x0);

  NumericVector res(x0.size());
  for(int i=0; i<res.size(); i++) res(i) = x0(i);
  //Rcpp::Rcout << "argmin      " << x0.transpose() << std::endl;
  //Rcpp::Rcout << "f in argmin " << f(x0) << std::endl;
  return res;
}

// [[Rcpp::export]]
NumericVector example_NewtonDescentSolver(NumericVector x, NumericVector mu, NumericMatrix sigma) {
  int k = mu.size();
  Eigen::VectorXd x_(x.size());
  Eigen::VectorXd mu_(mu.size());
  Eigen::MatrixXd sigma_(sigma.nrow(), sigma.ncol());
  for(int i=0; i<x.size(); i++) x_(i) = x(i);
  for(int i=0; i<mu.size(); i++) mu_(i) = mu(i);
  for(int i=0; i<sigma.nrow(); i++) for(int j=0; j<sigma.ncol(); j++) sigma_(i, j) = sigma(i, j);

  typedef ConditionalMode<double> TConditionalMode;
  typedef typename TConditionalMode::TVector TVector;
  typedef typename TConditionalMode::MatrixType MatrixType;

  TConditionalMode f(x_, mu_, sigma_);
  cppoptlib::NewtonDescentSolver<TConditionalMode> solver;
  Eigen::VectorXd x0(mu.size());
  for(int i=0; i<mu.size(); i++){
    x0(i) = log(x(i)+0.5) - log(x(k)+0.5);
  }
  solver.minimize(f, x0);

  NumericVector res(x0.size());
  for(int i=0; i<res.size(); i++) res(i) = x0(i);
  //Rcpp::Rcout << "argmin      " << x0.transpose() << std::endl;
  //Rcpp::Rcout << "f in argmin " << f(x0) << std::endl;
  return res;
}
