///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
#include <Eigen/SVD>

VectorXd lsqsvd(const MatrixXd &A, const VectorXd &b) {
  // Compute economical SVD, compare \cref{cpp:decompositions}
  Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  VectorXd sv = svd.singularValues();
  unsigned int r = svd.rank(); // Numerical rank, default tolerance
  MatrixXd U = svd.matrixU(), V = svd.matrixV();
  // $\cob{\Vx^{\dagger} = \VV_{1}\Sigmabf_{r}^{-1}\VU_{1}^{\herm}\Vb}$, see \eqref{lsq:svdsol}
  return V.leftCols(r) * (sv.head(r).cwiseInverse().asDiagonal() *
                          (U.leftCols(r).adjoint() * b));
}
/* SAM_LISTING_END_0 */

// Conversion into diagonal matrix could be replaced with componentwise scaling.

/* SAM_LISTING_BEGIN_1 */
VectorXd lsqsvd_eigen(const MatrixXd &A, const VectorXd &b) {
  Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  return svd.solve(b);
}
/* SAM_LISTING_END_1 */
