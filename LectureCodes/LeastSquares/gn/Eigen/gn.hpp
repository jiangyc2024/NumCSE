#include <Eigen/Dense>
#include <Eigen/QR>
using Eigen::VectorXd; 
using Eigen::MatrixXd;

template <class Function, class Jacobian>
VectorXd gn(const VectorXd& init, const Function& F, const Jacobian& J, const double tol) {
  VectorXd x = init;
  VectorXd s = J(x).householderQr().solve(F(x)); // \label{gn:2}
  x = x - s;
  while (s.norm() > tol * x.norm()) { // \label{gn:term}
    s = J(x).householderQr().solve(F(x)); // \label{gn:5}
    x = x - s;
  }

  return x;
}
