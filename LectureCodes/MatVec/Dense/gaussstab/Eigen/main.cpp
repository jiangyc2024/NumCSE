///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Dense>
#include <figure/figure.hpp>

using namespace std;
using namespace Eigen;

int main() {
  /* SAM_LISTING_BEGIN_0 */
  int n = 10;
  // Initialize vectors $\cob{\Vu}$ and $\cob{\Vv}$
  VectorXd u = VectorXd::LinSpaced(n, 1, n) / 3.0;
  VectorXd v = u.cwiseInverse().array() *
               VectorXd::LinSpaced(n, 1, n)
                   .unaryExpr([](double x) { return pow(-1, x); })
                   .array();
  VectorXd x = VectorXd::Ones(n);
  double nx = x.lpNorm<Infinity>();
  VectorXd expo = VectorXd::LinSpaced(19, -5, -14);
  Eigen::MatrixXd res(expo.size(), 4);
  for (int i = 0; i <= expo.size(); ++i) {
    // Build coefficient matrix $\cob{\VA}$
    double epsilon = std::pow(10, expo(i));
    MatrixXd A = u * v.transpose() + epsilon * MatrixXd::Identity(n, n);
    VectorXd b = A * x;               // right-hand-side vector
    double nb = b.lpNorm<Infinity>(); // maximum norm
    VectorXd xt = A.lu().solve(b);    // Gaussian elimination
    VectorXd r = b - A * xt;          // \com{residual vector}
    res(i, 0) = epsilon;
    res(i, 1) = (x - xt).lpNorm<Infinity>() / nx;
    res(i, 2) = r.lpNorm<Infinity>() / nb;
    // L-infinity condition number
    res(i, 3) = A.inverse().cwiseAbs().rowwise().sum().maxCoeff() *
                A.cwiseAbs().rowwise().sum().maxCoeff();
  }
  /* SAM_LISTING_END_0 */
  // Plotting
  mgl::Figure fig1;
  fig1.setlog(true, true);
  fig1.plot(res.col(0), res.col(1), " +r-").label("relative error");
  fig1.plot(res.col(0), res.col(2), " *b-").label("relative residual");
  fig1.xlabel("\\epsilon"); // not nice
  fig1.legend(0.05, 0.5);
  fig1.save("gausstab");
  return 0;
}
