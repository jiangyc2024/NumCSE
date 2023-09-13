#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>

namespace denseSolve {


using Eigen::MatrixXd;
using Eigen::VectorXd;


/* SAM_LISTING_BEGIN_0 */
// Gaussian elimination with \com{partial pivoting}, \cref{cpp:gepiv}
inline void lu_solve(const MatrixXd &A, const VectorXd &b, VectorXd &x) {
  x = A.lu().solve(b); // 'lu()' is short for 'partialPivLu()'
}

// Gaussian elimination with total pivoting
inline void fullpivlu_solve(const MatrixXd &A, const VectorXd &b, VectorXd &x) {
  x = A.fullPivLu().solve(b); // total pivoting
}

// An elimination solver based on Householder transformations
inline void qr_solve(const MatrixXd &A, const VectorXd &b, VectorXd &x) {
  const Eigen::HouseholderQR<MatrixXd> solver(A); // see \cref{sec:QR}
  x = solver.solve(b);
}

// Use singular value decomposition (SVD)
inline void svd_solve(const MatrixXd &A, const VectorXd &b, VectorXd &x) {
  // SVD based solvers, see \cref{sec:SVD}
  x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
}
/* SAM_LISTING_END_0 */


} //namespace denseSolve