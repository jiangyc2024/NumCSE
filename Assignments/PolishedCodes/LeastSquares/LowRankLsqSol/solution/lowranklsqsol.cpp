/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: December 2022
 */

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace LowRankLsq {
struct LowRankMatrix {
  LowRankMatrix(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B)
      : A_(A), B_(B) {
    if (A.cols() != B.cols()) {
      throw std::runtime_error("A,B size mismatch");
    }
  }
  operator Eigen::MatrixXd() const { return A_ * (B_.transpose()); }
  Eigen::MatrixXd A_;  // First matrix factor
  Eigen::MatrixXd B_;  // Transpose of second matrix factor
};

/* SAM_LISTING_BEGIN_0 */
unsigned int rank(const Eigen::MatrixXd &A,
                  double tol = std::numeric_limits<double>::epsilon()) {
  if (A.norm() == 0) return 0;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_object(A);
  const Eigen::VectorXd sv = svd_object.singularValues();
  const unsigned int n = sv.size();
  unsigned int r = 0;
  while ((r < n) && (sv(r) >= sv(0) * A.rows() * A.cols() * tol)) r++;
  return r;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
unsigned int rank(const LowRankMatrix &M,
                  double tol = std::numeric_limits<double>::epsilon()) {
  // Maximal rank of matrix M
  unsigned int p = M.A_.cols();
  if ((M.A_.norm() == 0.0) || (M.B_.norm() == 0.0)) return 0;
  // Compute singular values efficiently
  // QR decompositions of matrix factors
  Eigen::HouseholderQR<Eigen::MatrixXd> qrA(M.A_);
  Eigen::HouseholderQR<Eigen::MatrixXd> qrB(M.B_);
  // Extract R factors and form their product
  const Eigen::MatrixXd RA =
      qrA.matrixQR().block(0, 0, p, p).template triangularView<Eigen::Upper>();
  const Eigen::MatrixXd RB =
      qrB.matrixQR().block(0, 0, p, p).template triangularView<Eigen::Upper>();
  // Compute SVD of $\cob{\VR_A\VR_B^\top}$, singular values only
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_object(RA * (RB.transpose()));
  const Eigen::VectorXd sv = svd_object.singularValues();
  const unsigned int n = sv.size();
  unsigned int r = 0;
  while ((r < n) && (sv(r) >= sv(0) * M.A_.rows() * M.B_.rows() * tol)) r++;
  return r;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_7 */
Eigen::VectorXd lowRankLsqSol(
    const LowRankMatrix &M, const Eigen::VectorXd &b,
    double tol = std::numeric_limits<double>::epsilon()) {
  assert(M.B_.rows() == b.size());
  // Maximal rank of matrix M
  unsigned int p = M.A_.cols();
  // Economical (!) QR decompositions of matrix factors
  Eigen::HouseholderQR<Eigen::MatrixXd> qrA(M.A_);  // \Label[line]{lsv:qrA}
  Eigen::HouseholderQR<Eigen::MatrixXd> qrB(M.B_);  // \Label[line]{lsv:qrB}
  // Extract R factors and form their product
  const Eigen::MatrixXd RA =  // \Label[line]{lsv:RA}
      qrA.matrixQR().block(0, 0, p, p).template triangularView<Eigen::Upper>();
  const Eigen::MatrixXd RB =  // \Label[line]{lsv:RB}
      qrB.matrixQR().block(0, 0, p, p).template triangularView<Eigen::Upper>();
  // Compute SVD of $\cob{\VR_A\VR_B^\top}$
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_object(
      RA * (RB.transpose()),
      Eigen::ComputeThinU | Eigen::ComputeThinV);  // \Label[line]{lsv:svd}
  // Determine numerical rank
  const Eigen::VectorXd sv = svd_object.singularValues();
  const unsigned int n = sv.size();
  unsigned int r = 0;
  while ((r < n) && (sv(r) >= sv(0) * M.A_.rows() * M.B_.rows() * tol)) r++;
  // Compute SVD factors
  const Eigen::MatrixXd U =
      qrA.householderQ() *
    (Eigen::MatrixXd::Identity(M.A_.rows(), p) * svd_object.matrixU());  // \Label[line]{lsv:1}

  const Eigen::MatrixXd V =
      qrB.householderQ() *
      (Eigen::MatrixXd::Identity(M.B_.rows(), p) * svd_object.matrixV());  // \Label[line]{lsv:2}
  return V.leftCols(r) *
         (sv.head(r).cwiseInverse().asDiagonal() *  // \Label[line]{lsv:x}
          (U.leftCols(r).transpose() * b));
}
/* SAM_LISTING_END_7 */

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd> svd(
    const LowRankMatrix &M) {
  // Maximal rank of matrix M
  unsigned int p = M.A_.cols();
  // QR decompositions of matrix factors
  Eigen::HouseholderQR<Eigen::MatrixXd> qrA(M.A_);
  Eigen::HouseholderQR<Eigen::MatrixXd> qrB(M.B_);
  // Extract R factors and form their product
  const Eigen::MatrixXd RA =
      qrA.matrixQR().block(0, 0, p, p).template triangularView<Eigen::Upper>();
  const Eigen::MatrixXd RB =
      qrB.matrixQR().block(0, 0, p, p).template triangularView<Eigen::Upper>();
  // Compute SVD of $\cob{\VR_A\VR_B^\top}$
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_object(
      RA * (RB.transpose()), Eigen::ComputeThinU | Eigen::ComputeThinV);
  // Assemble SVD factors
  return {qrA.householderQ() * (Eigen::MatrixXd::Identity(M.A_.rows(), p) *
                                svd_object.matrixU()),
          svd_object.singularValues(),
          qrB.householderQ() * (Eigen::MatrixXd::Identity(M.B_.rows(), p) *
                                svd_object.matrixV())};
}

bool test_svd(const LowRankMatrix &M) {
  Eigen::MatrixXd M_dense = M;
  auto [U, s, V] = svd(M);
  const unsigned int p = U.cols();
  // Check sizes
  if ((p != V.cols()) || (p != s.size())) {
    std::cerr << "Size mismatch: U = " << U.rows() << "x" << U.cols()
              << ", V= " << V.rows() << "x" << V.cols() << ", #s = " << s.size()
              << std::endl;
    return false;
  }

  // Check orthonormality of columns of U and V
  if ((U.transpose() * U - Eigen::MatrixXd::Identity(p, p)).norm() > 1.0E-9) {
    std::cerr << "U does not have orthonormal columns" << std::endl;
    return false;
  }
  if ((V.transpose() * V - Eigen::MatrixXd::Identity(p, p)).norm() > 1.0E-9) {
    std::cerr << "V does not have orthonormal columns" << std::endl;
    return false;
  }
  // Check matrix factorization
  if ((U * s.asDiagonal() * (V.transpose()) - M_dense).norm() >
      1.0E-9 * M_dense.norm()) {
    std::cerr << "Factorization failure" << std::endl;
    return false;
  }
  return true;
}

}  // namespace LowRankLsq

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "NumCSE code: C++ wrapper class for least-squares solution of "
               "low-rank LSE"
            << std::endl;

  const int m = 10;
  const int n = 7;
  const int p = 2;
  // generate test matrices
  Eigen::MatrixXd A(m, p), B(n, p);
  for (int l = 0; l < m; ++l) {
    A(l, 0) = 1.0 + l;
    A(l, 1) = A(l, 0);  // 1.0 + 2 * l;
  }
  for (int l = 0; l < n; ++l) {
    B(l, 0) = 2.0 - l;
    B(l, 1) = 1.0 + 3 * l;
  }
  LowRankLsq::LowRankMatrix M(A, B);
  Eigen::MatrixXd M_dense{M};
  std::cout << "Test of low-rank SVD: "
            << (LowRankLsq::test_svd(M) ? "success" : "failure") << std::endl;
  Eigen::JacobiSVD<Eigen::MatrixXd> svdM(M_dense);
  std::cout << "SVD rank = " << svdM.rank() << " <-> AB^T-rank = " << rank(M)
            << std::endl;
  return 0;
}
