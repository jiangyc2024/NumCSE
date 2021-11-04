#ifndef POLYVAL_HPP
#define POLYVAL_HPP

#include <Eigen/Dense>

/* evaluate polynomial using horner scheme
 * coefficients in p must be given like:
 * f(x) = p0*x^n + p1*x^(n-1) + ... + pn
 *
 * IN:  p = coefficients as stated above
 *      x = Eigen::Eigen::VectorXd of points at which to evaluate
 * OUT: y = f(x)                                           */
void polyval(const Eigen::VectorXd& p, const Eigen::VectorXd& x,
             Eigen::VectorXd& y) {
  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
  y = p(0) * ones;

  for (unsigned i = 1; i < p.size(); ++i) {
    y = x.cwiseProduct(y) + p(i) * ones;
  }
}

/* does the same as polyval above but can be called as
 *   Eigen::VectorXd y = polyval(p, x);
 * instead of
 *   Eigen::VectorXd y;
 *   polyval(p, x, y);
 *
 * IN:  p = coefficients as stated above
 *      x = Eigen::Eigen::VectorXd of points at which to evaluate
 * OUT: y = f(x)                                           */
Eigen::VectorXd polyval(const Eigen::VectorXd& p, const Eigen::VectorXd& x) {
  Eigen::VectorXd y;
  polyval(p, x, y);
  return y;
}

#endif
