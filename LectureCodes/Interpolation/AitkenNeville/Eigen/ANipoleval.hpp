///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): J. Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

# include <Eigen/Dense>

using Eigen::VectorXd;
/* SAM_LISTING_BEGIN_0 */
// Aitken-Neville algorithm for evaluation of interpolating polynomial
// IN:  t, y: (vectors of) interpolation data points 
//      x: (single) evaluation point
// OUT: value of interpolant in x
double ANipoleval(const Eigen::VectorXd& t, Eigen::VectorXd y, const double x) {
  for (int i = 1; i < y.size(); ++i) {
    // Careful: int loop index required for comparison >= 0 !
    for (int k = i - 1; k >= 0; --k) {
      // Recursion \eqref{eq:ipolrec}
      y[k] = y[k + 1] + (y[k + 1] - y[k])*(x - t[i])/(t[i] - t[k]);
    }
  }
  return y[0];
}
/* SAM_LISTING_END_0 */
