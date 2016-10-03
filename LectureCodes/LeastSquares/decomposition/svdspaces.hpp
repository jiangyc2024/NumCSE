///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <limits>
#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;

# include <Eigen/SVD>
const double EPS = std::numeric_limits<double>::epsilon();

/* SAM_LISTING_BEGIN_0 */
// Computation of the numerical rank of a non-zero matrix by means of SVD
// A complete implementation would have to test for \Blue{$\VA\approx\Vzero$}
MatrixXd::Index rank_ncse(const MatrixXd &A, double tol = EPS) {
  Eigen::JacobiSVD<MatrixXd> svd(A);
  const VectorXd sv = svd.singularValues(); // Get \com{sorted} singular values as vector
  MatrixXd::Index r = 0;
  while (sv(r) >= sv(0)*tol) r++; // Test \com{relative} size of singular values
  return r;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
// Computation of the numerical rank of a matrix by means of SVD
MatrixXd::Index rank_eigen(const MatrixXd &A, double tol = EPS) {
  return A.jacobiSvd().setThreshold(tol).rank();
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
// Computation of an ONB of the kernel of a matrix
MatrixXd nullspace(const MatrixXd &A, double tol = EPS) {
  using index_t = MatrixXd::Index;
  Eigen::JacobiSVD<MatrixXd> svd(A,Eigen::ComputeFullV);
  index_t r = svd.setThreshold(tol).rank();
  // Rightmost columns of \Blue{$\VV$} provide ONB of \Blue{$\Kern(\VA)$}
  MatrixXd Z = svd.matrixV().rightCols(A.cols()-r); 
  return Z; 
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
// Computation of an ONB of the kernel of a matrix
MatrixXd rangespace(const MatrixXd &A, double tol = EPS) {
  using index_t = MatrixXd::Index;
  Eigen::JacobiSVD<MatrixXd> svd(A,Eigen::ComputeFullV);
  index_t r = svd.setThreshold(tol).rank();
  // r left columns of \Blue{$\VU$} provide ONB of \Blue{$\Range(\VA)$}
  return svd.matrixV().leftCols(r); 
}
/* SAM_LISTING_END_3 */




// Important! If replacing 'ComputeFull' by 'ComputeThin', 
// the methods matrixU() and matrixV() will return the 
// matrices U and V but not extended to an orthonormal basis!
