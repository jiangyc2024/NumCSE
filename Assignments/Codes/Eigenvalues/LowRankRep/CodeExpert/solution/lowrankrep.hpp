#ifndef LOWRANKREP_HPP
#define LOWRANKREP_HPP

////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <iostream>

using namespace Eigen;

/* @brief Factorize matrix $X$ into $X = AB'$
 * @param[in] X An $m \times n$ matrix of rank at most $k$
 * @param[in] k Rank of matrix $X$
 * @param[out] A The $m \times k$ matrix from the decomposition of $X$
 * @param[out] B The $n \times k$ matrix from the decomposition of $X$
 */
/* SAM_LISTING_BEGIN_0 */
std::pair<MatrixXd, MatrixXd> factorize_X_AB(const MatrixXd &X,
                                             const unsigned int k) {
  size_t m = X.rows(), n = X.cols();

  double tol = 1e-6; // Tolerance for numerical rank: \lref{ex:svdrank}
  // TO DO (4-8.b) Define A,B with k columns such that X = A*B^T
  // Hint: warnings can be displayed with std::cerr, error messages with assert
  // START
  assert(k <= std::min(m, n) && "Rank k cannot be larger than dimensions of X");

  JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
  // With Eigen::svd you can ask for thin $\VU$ or $\VV$ to be computed.
  // In case of a rectangular $m \times n$ matrix,
  // with $j$ the smaller value among $m$ and $n$,
  // there can only be at most $j$ singular values.
  // The remaining columns of $\VU$ and $\VV$ do not correspond
  // to actual singular vectors and are not computed in thin format.

  VectorXd s = svd.singularValues();

  if (k + 1 <= std::min(m, n) && s(k) > tol * s(0)) {
    // 1.\ condition checks if there is a singular value after $k$.
    // ($\VS$ cannot be larger than the dimensions of $\VX$.)
    // 2.\ condition checks if such singular value
    // (by definition non-negative) is "numerically" different from 0.
    std::cerr << "Rank of matrix X is greater than required rank k"
              << std::endl;
  }
  // We know that k <= std::min(m,n) by the assert
  MatrixXd A = svd.matrixU().leftCols(k) * s.asDiagonal();
  MatrixXd B = svd.matrixV().leftCols(k);
  // END
  return std::make_pair(A, B);
}
/* SAM_LISTING_END_0 */

/* @brief Compute the SVD of $AB'$
 * @param[in] A An $m \times k$ matrix
 * @param[in] B An $n \times k$ matrix
 * @param[out] U The $n \times k$ matrix from the SVD of $AB'$
 * @param[out] S The $k \times k$ diagonal matrix of sing. vals of $AB'$
 * @param[out] V The $n \times k$ matrix from the SVD of $AB'$
 */
/* SAM_LISTING_BEGIN_1 */
std::tuple<MatrixXd, MatrixXd, MatrixXd> svd_AB(const MatrixXd &A,
                                                const MatrixXd &B) {

  assert(A.cols() == B.cols() &&
         "Matrices A and B should have the same column number");
  size_t m = A.rows();
  size_t n = B.rows();
  size_t k = A.cols();
  // TO DO (4-8.d): Define the tuple of
  // matrices U,S,V of the svd decomposition of A*B^T
  // Hint: You can define a tuple of objects using std::make_tuple
  // START
  // QA: m x k; QB: n x k; RA,RB k x k
  HouseholderQR<MatrixXd> QRA = A.householderQr();
  MatrixXd QA = QRA.householderQ() * MatrixXd::Identity(m, std::min(m, k));
  MatrixXd RA = MatrixXd::Identity(std::min(m, k), m) *
                QRA.matrixQR().triangularView<Upper>();
  HouseholderQR<MatrixXd> QRB = B.householderQr();
  MatrixXd QB = QRB.householderQ() * MatrixXd::Identity(n, std::min(n, k));
  MatrixXd RB = MatrixXd::Identity(std::min(n, k), n) *
                QRB.matrixQR().triangularView<Upper>();

  // U,V: k x k
  JacobiSVD<MatrixXd> svd(RA * RB.transpose(), ComputeFullU | ComputeFullV);
  // Thin matrices are unnecessary here as $RA*RB'$ is
  // a square $k \times k$ matrix!
  // Moreover, if you had computed the thin matrices $U$ and $V$
  // of the direct SVD of $AB'$,
  // you would still have dealt with the minimum between $m$ and $n$,
  // without exploiting $k << m, n$!

  VectorXd s = svd.singularValues();
  MatrixXd U = svd.matrixU();
  MatrixXd V = svd.matrixV();
  // U: m x k; V: n x k
  U = QA * U;
  V = QB * V;

  return std::make_tuple(U, s.asDiagonal(), V);
  // END
}
/* SAM_LISTING_END_1 */

/* @brief Find $Az$ and $Bz$ such that
 * $Z = Az*Bz'$ is the best approximation of $AxBx'+AyBy'$ with rank $k$
 * @param[in] Ax An $m \times k$ matrix
 * @param[in] Ay An $m \times k$ matrix
 * @param[in] Bx An $n \times k$ matrix
 * @param[in] By An $n \times k$ matrix
 * @param[out] Az The $m \times k$ matrix to form $Az*Bz'$
 * @param[out] Bz The $n \times k$ matrix to form $Az*Bz'$
 */
/* SAM_LISTING_BEGIN_2 */
std::pair<MatrixXd, MatrixXd> rank_k_approx(const MatrixXd &Ax,
                                            const MatrixXd &Ay,
                                            const MatrixXd &Bx,
                                            const MatrixXd &By) {

  assert(Ax.rows() == Ay.rows() && Ax.cols() == Ay.cols() &&
         "Matrices Ax and Ay should have the same dimensions");
  assert(Bx.rows() == By.rows() && Bx.cols() == By.cols() &&
         "Matrices Bx and By should have the same dimensions");
  // TO DO (4-8.h): Use the function svd_AB with the appropriate A,B and
  // return the pair Az, Bz of matrices of rank k with Z = Az*Bz
  // Hint: you can define pairs of objects using std::make_pair
  // Hint: use std::get<i>(my_tuple) to extract the i-th element of a tuple
  // START
  MatrixXd A(Ax.rows(), Ax.cols() + Ay.cols());
  A << Ax, Ay;
  MatrixXd B(Bx.rows(), Bx.cols() + By.cols());
  B << Bx, By;
  // U: m x 2k; S: 2k x 2k; V: n x 2k
  MatrixXd U, S, V;
  std::tuple<MatrixXd, MatrixXd, MatrixXd> svd_tuple = svd_AB(A, B);

  size_t k = Ax.cols();
  MatrixXd Az = std::get<0>(svd_tuple).leftCols(k) * std::get<1>(svd_tuple);
  MatrixXd Bz = std::get<2>(svd_tuple).leftCols(k);

  return std::make_pair(Az, Bz);
  // END
}
/* SAM_LISTING_END_2 */

#endif
