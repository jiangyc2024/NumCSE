# include <Eigen/Dense>

using Eigen::VectorXd;

/* evaluation of a polynomial in monomial representation using Horner scheme
 * IN: p = vector of coefficients, dim = 1x(deg + 1)
 *     x = evaluation points 
 * OUT: y = polynomial evaluated at x */
void horner (const VectorXd& p, const VectorXd& x, VectorXd& y) {
  const VectorXd ones = VectorXd::Ones(x.size());
  y = p(0)*ones;

  for (unsigned i = 1; i < p.size(); ++i) {
    y = x.cwiseProduct(y) + p(i)*ones;
  }
}
