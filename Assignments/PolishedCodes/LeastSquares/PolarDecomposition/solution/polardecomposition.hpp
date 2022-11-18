#ifndef POLARDECOMPOSITIONHPP
#define POLARDECOMPOSITIONHPP

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <iostream>

class PolarDecomposition {
 public:
  explicit PolarDecomposition(const Eigen::MatrixXd &X) { initialize(X); }
  PolarDecomposition(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);
  PolarDecomposition(const PolarDecomposition &) = default;
  ~PolarDecomposition() = default;

  // Left multiplication of M with the Q-factor of the polar decomposition
  void applyQ(Eigen::MatrixXd &Y) const { Y.applyOnTheLeft(Q_); }
  // Left multiplication of M with the M-factor of the polar decomposition
  void applyM(Eigen::MatrixXd &Y) const { Y.applyOnTheLeft(M_); }

  int Qcols() { return Q_.cols(); }
  int Mcols() { return M_.cols(); }

 private:
  void initialize(const Eigen::MatrixXd &X);
  Eigen::MatrixXd Q_;  // factor Q
  Eigen::MatrixXd M_;  // factor M
};

/* SAM_LISTING_BEGIN_1 */
void PolarDecomposition::initialize(const Eigen::MatrixXd &X) {
  assert(X.rows() >= X.cols());
  // TODO: (3-12.c) Implement a method to initialize the data members Q_ and M_
  // corresponding to Q and M in Theorem 0.3.1, where X = QM
  // START
  // Compute SVD
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      X, Eigen::ComputeThinU | Eigen::ComputeFullV);
  // Extract factor matrices
  const Eigen::MatrixXd &V{svd.matrixV()};
  const Eigen::MatrixXd &U{svd.matrixU()};
  const auto Sigma{svd.singularValues().asDiagonal()};
  // Form factors of polar decomposition
  Q_ = U * V.transpose();
  M_ = V * Sigma * V.transpose();
  // END
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_7 */
PolarDecomposition::PolarDecomposition(const Eigen::MatrixXd &A,
                                       const Eigen::MatrixXd &B) {
  const long m = A.rows();  // No. of rows of $\cob{\VX}$
  const long n = B.rows();  // No. of columns of $\cob{\VX}$
  const long k = A.cols();  // Maximal rank of $\cob{\VX}$
  // We assume $\cob{k \leq n \leq m}$
  assert(m >= n);
  assert(k < n);
  assert(B.cols() == k);
  // TODO: (3-12.d) Implement a method to initialize the data members Q_ and M_
  // for X := AB^T = QM, with optimal complexity
  // START
  // Compute QR-decompositions in an encoded way, see \lref{par:ecovsfull}
  Eigen::HouseholderQR<Eigen::MatrixXd> QRA(A);  // cost = $\cob{O(mk^2)}$
  Eigen::HouseholderQR<Eigen::MatrixXd> QRB(B);  // cost = $\cob{O(nk^2)}$
  const Eigen::MatrixXd RA{
      QRA.matrixQR().block(0, 0, k, k).template triangularView<Eigen::Upper>()};
  const Eigen::MatrixXd RB{
      QRB.matrixQR().block(0, 0, k, k).template triangularView<Eigen::Upper>()};
  // SVD of small $\cob{k\times k}$-matrix $\cob{\VR_A\VR_B^{\top}}$
  // cost = $\cob{O(k^3)}$
  Eigen::JacobiSVD<Eigen::MatrixXd> svdh(
      RA * RB.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
  // Extract kxk orthogonal factor matrices
  const Eigen::MatrixXd &Vh{svdh.matrixV()};
  const Eigen::MatrixXd &Uh{svdh.matrixU()};
  // Build auxiliary mxn block matrix
  Eigen::MatrixXd W{Eigen::MatrixXd::Zero(m, n)};  // cost = $\cob{O(mn)}$
  W.block(0, 0, k, k) = Uh * Vh.transpose();       // cost = $\cob{O(k^3)}$
  W.block(k, k, n - k, n - k) = Eigen::MatrixXd::Identity(n - k, n - k);
  // Compute the Q-factor, cost = $\cob{O(kmn+kn^2)}$
  Q_ =
      (QRB.householderQ() * ((QRA.householderQ() * W).transpose())).transpose();
  // Small kxk matrix containing singular values on diagonal
  const auto Sigma{svdh.singularValues().asDiagonal()};
  // Form M-factor of polar decomposition
  // Auxiliary matrices
  Eigen::MatrixXd S{Eigen::MatrixXd::Zero(n, k)};
  S.block(0, 0, k, k) = Vh * Sigma * Vh.transpose();  // cost = $\cob{O(k^3)}$
  Eigen::MatrixXd T{Eigen::MatrixXd::Zero(n, n)};     // cost = $\cob{O(n^2)}$
  T.block(0, 0, k, n) =
      (QRB.householderQ() * S).transpose();   // cost = $\cob{O(kn^2)}$
  M_ = (QRB.householderQ() * T).transpose();  // cost = $\cob{O(kn^2)}$
  // END
}
/* SAM_LISTING_END_7 */

#endif
