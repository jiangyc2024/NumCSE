#ifndef LOWRANKREP_HPP
#define LOWRANKREP_HPP

////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////

#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <utility>

/**
 * @brief Factorize matrix $X$ into $X = AB'$
 *
 * @param X An $m \times n$ matrix of rank at most $k$
 * @param k Rank of matrix $X$
 * @return std::pair<Eigen::MatrixXd, Eigen::MatrixXd> pair (A, B) with A: The
 * $m \times k$ matrix from the decomposition of $X$ and B: The $n \times k$
 * matrix from the decomposition of $X$
 */
/* SAM_LISTING_BEGIN_0 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> factorize_X_AB(
    const Eigen::MatrixXd& X, const unsigned int k) {
  std::size_t m = X.rows(), n = X.cols();
  Eigen::MatrixXd A, B;
  constexpr double tol =
      1e-6;  // Tolerance for numerical rank: \lref{ex:svdrank}
  // TODO: (3-8.b) Define A,B with k columns such that X = A*B^T
  // Hint: warnings can be displayed with std::cerr, error messages with assert
  // START
  assert(k <= std::min(m, n) && "Rank k cannot be larger than dimensions of X");

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      X, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // With Eigen::svd you can ask for thin $\VU$ or $\VV$ to be computed.
  // In case of a rectangular $m \times n$ matrix,
  // with $j$ the smaller value among $m$ and $n$,
  // there can only be at most $j$ singular values.
  // The remaining columns of $\VU$ and $\VV$ do not correspond
  // to actual singular vectors and are not computed in thin format.

  Eigen::VectorXd s = svd.singularValues();

  if (k + 1 <= std::min(m, n) && s(k) > tol * s(0)) {
    // 1.\ condition checks if there is a singular value after $k$.
    // ($\VS$ cannot be larger than the dimensions of $\VX$.)
    // 2.\ condition checks if such singular value
    // (by definition non-negative) is "numerically" different from 0.
    std::cerr << "Rank of matrix X is greater than required rank k"
              << std::endl;
  }
  // We know that k <= std::min(m,n) by the assert
  A = svd.matrixU().leftCols(k) * s.asDiagonal();
  B = svd.matrixV().leftCols(k);
  // END
  return std::make_pair(A, B);
}
/* SAM_LISTING_END_0 */

/**
 * @brief Compute the SVD of $AB'$
 *
 * @param A An $m \times k$ matrix
 * @param B An $n \times k$ matrix
 * @return std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> tuple
 * (U, S, V) with U: The $n \times k$ matrix from the SVD of $AB'$, S: The $k
 * \times k$ diagonal matrix of sing. vals of $AB'$ and V: The $n \times k$
 * matrix from the SVD of $AB'$
 */
/* SAM_LISTING_BEGIN_1 */
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> svd_AB(
    const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
  assert(A.cols() == B.cols() &&
         "Matrices A and B should have the same column number");
  assert(A.cols() < A.rows() && A.cols() < B.rows() &&
         "Matrices A and B should have less columns than rows");
  std::size_t m = A.rows();
  std::size_t n = B.rows();
  std::size_t k = A.cols();
  Eigen::MatrixXd U, S, V;
  // TODO: (3-8.d) Define the tuple of
  // matrices U,S,V of the svd decomposition of A*B^T
  // Hint: You can define a tuple of objects using std::make_tuple
  // START
  // QA: m x k; QB: n x k; RA,RB k x k
  Eigen::HouseholderQR<Eigen::MatrixXd> QRA = A.householderQr();
  Eigen::MatrixXd QA =
      QRA.householderQ() * Eigen::MatrixXd::Identity(m, std::min(m, k));
  Eigen::MatrixXd RA = Eigen::MatrixXd::Identity(std::min(m, k), m) *
                       QRA.matrixQR().triangularView<Eigen::Upper>();
  Eigen::HouseholderQR<Eigen::MatrixXd> QRB = B.householderQr();
  Eigen::MatrixXd QB =
      QRB.householderQ() * Eigen::MatrixXd::Identity(n, std::min(n, k));
  Eigen::MatrixXd RB = Eigen::MatrixXd::Identity(std::min(n, k), n) *
                       QRB.matrixQR().triangularView<Eigen::Upper>();

  // U,V: k x k
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      RA * RB.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
  // Thin matrices are unnecessary here as $RA*RB'$ is
  // a square $k \times k$ matrix!
  // Moreover, if you had computed the thin matrices $U$ and $V$
  // of the direct SVD of $AB'$,
  // you would still have dealt with the minimum between $m$ and $n$,
  // without exploiting $k << m, n$!

  S = svd.singularValues().asDiagonal();
  U = svd.matrixU();
  V = svd.matrixV();
  // U: m x k; V: n x k
  U = QA * U;
  V = QB * V;
  // END
  return std::make_tuple(U, S, V);
}
/* SAM_LISTING_END_1 */

/**
 * @brief Find $Az$ and $Bz$ such that
 * $Z = Az*Bz'$ is the best approximation of $AxBx'+AyBy'$ with rank $k$
 *
 * @param Ax An $m \times k$ matrix
 * @param Ay An $m \times k$ matrix
 * @param Bx An $n \times k$ matrix
 * @param By An $n \times k$ matrix
 * @return std::pair<Eigen::MatrixXd, Eigen::MatrixXd> pair (Az, Bz) with Az:
 * The $m \times k$ matrix to form $Az*Bz'$ and Bz: The $n \times k$ matrix to
 * form $Az*Bz'$
 */
/* SAM_LISTING_BEGIN_2 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> rank_k_approx(
    const Eigen::MatrixXd& Ax, const Eigen::MatrixXd& Ay,
    const Eigen::MatrixXd& Bx, const Eigen::MatrixXd& By) {
  assert(Ax.rows() == Ay.rows() && Ax.cols() == Ay.cols() &&
         "Matrices Ax and Ay should have the same dimensions");
  assert(Bx.rows() == By.rows() && Bx.cols() == By.cols() &&
         "Matrices Bx and By should have the same dimensions");
  assert(Ax.cols() == Bx.cols() &&
         "Matrices Ax, Ay, Bx and By should have the same amount of columns");
  assert(Ax.cols() * 2 < Ax.rows() && Bx.cols() * 2 < Bx.rows() &&
         "2*k should not be bigger equal min(m, n)");

  Eigen::MatrixXd Az, Bz;
  // TODO: (3-8.h) Use the function svd_AB with the appropriate A,B and
  // return the pair Az, Bz of matrices of rank k with Z = Az*Bz
  // Hint: you can define pairs of objects using std::make_pair
  // Hint: use std::get<i>(my_tuple) to extract the i-th element of a tuple
  // START
  Eigen::MatrixXd A(Ax.rows(), Ax.cols() + Ay.cols());
  A << Ax, Ay;
  Eigen::MatrixXd B(Bx.rows(), Bx.cols() + By.cols());
  B << Bx, By;
  // U: m x 2k; S: 2k x 2k; V: n x 2k
  auto [U, S, V] = svd_AB(A, B);
  std::size_t k = Ax.cols();
  Az = U.leftCols(k) * S.topLeftCorner(k, k);
  Bz = V.leftCols(k);
  // END
  return std::make_pair(Az, Bz);
}
/* SAM_LISTING_END_2 */

#endif
