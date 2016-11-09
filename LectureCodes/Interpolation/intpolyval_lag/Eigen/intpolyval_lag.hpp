///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): J. Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <Eigen/Dense>
# include "lagrangepoly.hpp"

using Eigen::VectorXd;
/* SAM_LISTING_BEGIN_0 */
// Evaluation of the interpolation polynomials with Lagrange polynomials
// IN:  t = vector of interpolation nodes
//      y = values in t
//      x = evaluation points
// OUT  p will be used to save the values of the interpolant evaluated in x 
void intpolyval_lag(const VectorXd& t, const VectorXd& y, const VectorXd& x, VectorXd& p) {
  p = VectorXd::Zero(x.size());
  for (unsigned k = 0; k < t.size(); ++k) {
    // Compute values of k-th Lagrange polynomial in evaluation points
    VectorXd L; lagrangepoly(x, k, t, L);
    p += y(k)*L;
  }
}
/* SAM_LISTING_END_0 */
