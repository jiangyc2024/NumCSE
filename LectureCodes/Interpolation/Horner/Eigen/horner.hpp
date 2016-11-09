///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <Eigen/Dense>

using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
/* Efficient evaluation of a polynomial in monomial representation 
 * using the Horner scheme \eqref{intp:Horner}
 * IN: p = vector of \com{monomial coefficients}, length = degree + 1
 *     x = vector of evaluation points 
 * OUT: y = polynomial evaluated at x */
void horner(const VectorXd& p, const VectorXd& x, VectorXd& y) {
  const VectorXd::Index n = x.size();
  y.resize(n); y = p(0)*VectorXd::Ones(n);
  for (unsigned i = 1; i < p.size(); ++i) {
    y = x.cwiseProduct(y) + p(i)*VectorXd::Ones(n);
  }
}
/* SAM_LISTING_END_0 */
