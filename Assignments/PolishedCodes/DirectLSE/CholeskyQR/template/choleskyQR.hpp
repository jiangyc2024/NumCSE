#ifndef CHOLESKYQR_HPP
#define CHOLESKYQR_HPP

#include <Eigen/Dense>
#include <algorithm>

/**
 * @brief QR decomposition from Cholesky decomposition
 *
 * @param A An $m \times n$ matrix
 * @param R The upper triangular matrix from the QR decomposition of $A$
 * @param Q The orthogonal matrix from the QR decomposition of $A$
 */
/* SAM_LISTING_BEGIN_0 */
void CholeskyQR(const Eigen::MatrixXd& A, Eigen::MatrixXd& R,
                Eigen::MatrixXd& Q) {
  Eigen::MatrixXd AtA = A.transpose() * A;
  Eigen::LLT<Eigen::MatrixXd> L = AtA.llt();
  R = L.matrixL().transpose();
  Q = R.transpose()
          .triangularView<Eigen::Lower>()
          .solve(A.transpose())
          .transpose();
  // .triangularView() template member only accesses the triangular part
  // of a dense matrix and allows to easily solve linear problem
}
/* SAM_LISTING_END_0 */

/**
 * @brief Direct QR decomposition
 *
 * @param A An $m \times n$ matrix
 * @param R The upper triangular matrix from the QR decomposition of $A$
 * @param Q The orthogonal matrix from the QR decomposition of $A$
 */
/* SAM_LISTING_BEGIN_1 */
void DirectQR(const Eigen::MatrixXd& A, Eigen::MatrixXd& R,
              Eigen::MatrixXd& Q) {
  const size_t m = A.rows();
  const size_t n = A.cols();

  Eigen::HouseholderQR<Eigen::MatrixXd> QR = A.householderQr();
  Q = QR.householderQ() * Eigen::MatrixXd::Identity(m, std::min(m, n));
  R = Eigen::MatrixXd::Identity(std::min(m, n), m) *
      QR.matrixQR().triangularView<Eigen::Upper>();
  // If A: m x n, then Q: m x m and R: m x n.
  // If m > n, however, the extra columns of Q and extra rows of R are not
  // needed. Matlab returns this "economy-size" format calling "qr(A,0)", which
  // does not compute these extra entries. With the code above, Eigen is smart
  // enough to not compute the discarded vectors.
}
/* SAM_LISTING_END_1 */

#endif