# ifndef POLYVAL_HPP
# define POLYVAL_HPP

# include <Eigen/Dense>

using Eigen::VectorXd;

/* evaluate polynomial using horner scheme 
 * coefficients in p must be given like: 
 * f(x) = p0*x^n + p1*x^(n-1) + ... + pn
 *
 * IN:  p = coefficients as stated above
 *      x = Eigen::VectorXd of points at which to evaluate
 * OUT: y = f(x)                                           */
void polyval(const VectorXd& p, const VectorXd& x, VectorXd& y) {
  const VectorXd ones = VectorXd::Ones(x.size());
  y = p(0)*ones;

  for (unsigned i = 1; i < p.size(); ++i) {
    y = x.cwiseProduct(y) + p(i)*ones;
  }
}

/* does the same as polyval above but can be called as
 *   VectorXd y = polyval(p, x);
 * instead of
 *   VectorXd y;
 *   polyval(p, x, y);
 *
 * IN:  p = coefficients as stated above
 *      x = Eigen::VectorXd of points at which to evaluate
 * OUT: y = f(x)                                           */
VectorXd polyval(const VectorXd& p, const VectorXd& x) {
  VectorXd y;
  polyval(p, x, y);
  return y;
}

# endif
