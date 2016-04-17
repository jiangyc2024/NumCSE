# include <Eigen/Dense>
# include "lagrangepoly.hpp"

using Eigen::VectorXd;
/* Evaluation of the interpolation polynomials with Lagrange polynomials
 * IN:  t = nodes
 *      y = values in t
 *      x = evaluation points
 * OUT: p = interpolant evaluated in x   */
void intpolyval_lag(const VectorXd& t, const VectorXd& y, const VectorXd& x, VectorXd& p) {
  p = VectorXd::Zero(x.size());
  for (unsigned k = 0; k < t.size(); ++k) {
    VectorXd L;
    lagrangepoly(x, k, t, L);
    p += y(k)*L;
  }
}
