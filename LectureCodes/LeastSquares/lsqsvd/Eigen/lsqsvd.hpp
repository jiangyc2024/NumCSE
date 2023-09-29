///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace lsqsvd {

using Eigen::MatrixXd;
using Eigen::VectorXd;

inline
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd lsqsvd(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
  // Compute economical SVD, compare \cref{cpp:decompositions}
  const Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  const Eigen::VectorXd & sv = svd.singularValues();
  const unsigned int r = svd.rank();  // Numerical rank, default tolerance
  // $\cob{\Vx^{\dagger} = \VV_{1}\Sigmabf_{r}^{-1}\VU_{1}^{\herm}\Vb}$, see \eqref{lsq:svdsol}
  const MatrixXd & U = svd.matrixU();
  const MatrixXd & V = svd.matrixV();
  
  return V.leftCols(r) * (sv.head(r).cwiseInverse().asDiagonal() *
                          (U.leftCols(r).adjoint() * b));
}
/* SAM_LISTING_END_0 */

// Conversion into diagonal matrix could be replaced with componentwise scaling.
// However, the expression template mechanism of Eigen should do this automatically.

inline
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd lsqsvd_eigen(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
  const Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  return svd.solve(b);
}
/* SAM_LISTING_END_1 */

inline
/**
 * @brief low-rank compression of a dense matrix via SVD
 * @param A real matrix of arbitrary size
 * @param k desired rank of compressed matrix
 * @return best rank-k compression of A
 * This function computes the best rank-k approximation of a given dense matrix
 * by truncating its economical SVD.
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd lowrankbestapprox(const Eigen::MatrixXd &A, unsigned int k) {
  // Compute economical SVD, compare \cref{cpp:decompositions}
  const Eigen::JacobiSVD<MatrixXd> svd(
      A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // Form matrix product $\cob{\VU_k\Sigmabf_k\VV_k}$.
  // Extract $\cob{\Sigmabf_k}$ as diagonal matrix of largest $k$ singular
  // values. \eigen provides singular values in decreasing order!
  return (svd.matrixU().leftCols(k)) *
         (svd.singularValues().head(k).asDiagonal()) *
         (svd.matrixV().leftCols(k).transpose());
}
/* SAM_LISTING_END_3 */


} // namespace lsqsvd