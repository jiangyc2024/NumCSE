///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>

using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// Efficient evaluation of a polynomial in monomial representation
// using the Horner scheme \eqref{intp:Horner}
// IN: p = vector of \com{monomial coefficients}, length = degree + 1
// (leading coefficient in p(0), \matlab convention \cref{rem:polyMatlab})
//     t = vector of evaluation points \Blue{$t_i$}
// OUT: vector of values: polynomial evaluated at \Blue{$t_i$}
Eigen::VectorXd horner(const Eigen::VectorXd &p, const Eigen::VectorXd &t) {
  const VectorXd::Index n = t.size();
  Eigen::VectorXd y{p[0] * VectorXd::Ones(n)};
  for (unsigned i = 1; i < p.size(); ++i)
    y = t.cwiseProduct(y) + p[i] * VectorXd::Ones(n);
  return y;
}
/* SAM_LISTING_END_0 */
