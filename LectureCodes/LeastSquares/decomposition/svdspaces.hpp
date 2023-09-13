///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <limits>
const double EPS = std::numeric_limits<double>::epsilon();

inline
/* SAM_LISTING_BEGIN_0 */
// Computation of the numerical rank of a non-zero matrix by means of
// singular value decomposition, cf. \eqref{eq:numrank}.
Eigen::Index rank_by_svd(const Eigen::MatrixXd &A, double tol = EPS) {
  if (A.norm() == 0) {
    return static_cast<Eigen::Index>(0);
  }
  const Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
  const Eigen::VectorXd & sv =
      svd.singularValues();  // Get \com{sorted} singular values as vector
  const Eigen::Index n = sv.size();
  Eigen::Index r = 0;
  // Test \com{relative} size of singular values
  while ((r < n) && (sv(r) >= sv(0) * tol)) {
    r++;
  }
  return r;
}
/* SAM_LISTING_END_0 */

inline
/* SAM_LISTING_BEGIN_1 */
// Computation of the numerical rank of a matrix by means of SVD
Eigen::Index rank_eigen(const Eigen::MatrixXd &A, double tol = EPS) {
  return A.jacobiSvd().setThreshold(tol).rank();
}
/* SAM_LISTING_END_1 */

inline
/* SAM_LISTING_BEGIN_2 */
// Computation of an ONB of the kernel of a matrix
Eigen::MatrixXd nullspace(const Eigen::MatrixXd &A, double tol = EPS) {
  using index_t = Eigen::Index;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
  const index_t r = svd.setThreshold(tol).rank();
  // Rightmost columns of \Blue{$\VV$} provide ONB of \Blue{$\Kern(\VA)$}
  Eigen::MatrixXd Z = svd.matrixV().rightCols(A.cols() - r);
  return Z;
}
/* SAM_LISTING_END_2 */

inline
/* SAM_LISTING_BEGIN_3 */
// Computation of an ONB of the image space of a matrix
Eigen::MatrixXd rangespace(const Eigen::MatrixXd &A, double tol = EPS) {
  using index_t = Eigen::Index;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU);
  const index_t r = svd.setThreshold(tol).rank();
  // r left columns of \Blue{$\VU$} provide ONB of \Blue{$\Range(\VA)$}
  return svd.matrixU().leftCols(r);
}
/* SAM_LISTING_END_3 */

// Important! If replacing 'ComputeFull' by 'ComputeThin',
// the methods matrixU() and matrixV() will return the
// matrices U and V but not extended to an orthonormal basis!
