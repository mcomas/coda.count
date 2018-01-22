#include <Rcpp.h>
#include "Eigen/Dense"
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include "cppoptlib/solver/bfgssolver.h"

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


using namespace cppoptlib;
using Eigen::VectorXd;

class Rosenbrock : public Problem<double> {
public:
  using typename cppoptlib::Problem<double>::Scalar;
  using typename cppoptlib::Problem<double>::TVector;
  void init() {

  }
  double value(const TVector &x) {
    const double t1 = (1 - x[0]);
    const double t2 = (x[1] - x[0] * x[0]);
    return   t1 * t1 + 100 * t2 * t2;
  }
  void gradient(const TVector &x, TVector &grad) {
    grad[0]  = -2 * (1 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2 * x[0]);
    grad[1]  = 200 * (x[1] - x[0] * x[0]);
  }
};

// [[Rcpp::export]]
VectorXd example02() {
  Rosenbrock f;
  BfgsSolver<Rosenbrock> solver;
  VectorXd x(2); x << -1, 2;
  solver.minimize(f, x);
  Rcpp::Rcout << "argmin      " << x.transpose() << std::endl;
  Rcpp::Rcout << "f in argmin " << f(x) << std::endl;
  return x;
}
