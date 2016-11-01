///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;

/* SAM_LISTING_BEGIN_0 */
// Initialization of a Vandermonde matrix \eqref{intp:VandermondeM}
// from interpolation points \Blue{$t_i$}.
MatrixXd vander(const VectorXd &t) {
  const VectorXd::Index n = t.size();
  MatrixXd V(n,n); V.col(0) = VectorXd::Ones(n); V.col(1) = t;
  // Store componentwise integer powers of point coordinate vector
  // into the columns of the Vandermonde matrix
  for (int j=2;j<n;j++) V.col(j) = (t.array().pow(j)).matrix();
  return V;
}
/* SAM_LISTING_END_0 */
