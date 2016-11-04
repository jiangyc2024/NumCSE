///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

# ifndef POLYFIT_HPP
# define POLYFIT_HPP

# include <Eigen/Dense>
# include <Eigen/QR>

using Eigen::VectorXd;
using Eigen::MatrixXd;

/* SAM_LISTING_BEGIN_0 */
// Solver for polynomial linear least squares data fitting problem
// data points passed in t and y, 'order' = degree + 1 
VectorXd polyfit(const VectorXd& t, const VectorXd& y, const unsigned& order) {
  // Initialize the coefficient matrix of \eqref{eq:polyfitlse}
  Eigen::MatrixXd A = Eigen::MatrixXd::Ones(t.size(),order + 1);
  for (unsigned j = 1; j < order + 1; ++j) 
    A.col(j) = A.col(j - 1).cwiseProduct(t);
  // Use \eigen's built-in least squares solver, see \cref{cpp:lsqsolveeigen}
  Eigen::VectorXd coeffs = A.householderQr().solve(y);
  // leading coefficients have low indices.
  return coeffs.reverse();
}
/* SAM_LISTING_END_0 */

# endif
