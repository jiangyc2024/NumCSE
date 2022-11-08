#ifndef CONCATMATLRC_HPP
#define CONCATMATLRC_HPP

/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#define _USE_MATH_DEFINES

/**
 * @brief computes the factors of the economical SVD of a matrix created by
 * concatenating two equal-size matrices given in low-rank factorization
 *
 * @param A1 matrix 1, left factor of low-rank factorization
 * @param B1 matrix 1, right factor of low-rank factorization
 * @param A2 matrix 2, left factor of low-rank factorization
 * @param B2 matrix 2, right factor of low-rank factorization
 * @return std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd> factors
 * of SVD
 */
/* SAM_LISTING_BEGIN_1 */
std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd> eco_svd_concat(
    const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
    const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2) {
  // Matrix dimensions and consistency checks
  const unsigned int m = A1.rows();
  assert(A2.rows() == m);
  const unsigned int n = B1.rows();
  assert(B2.rows() == n);
  const unsigned int k = A1.cols();
  assert((k == A2.cols()) && (k == B1.cols()) && (k == B2.cols()));
  // Matrices/vector to be returned by the function
  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(m, 2 * k);
  Eigen::VectorXd s = Eigen::VectorXd::Zero(2 * k);
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(2 * n, 2 * k);

  // TODO: (3-15.b) Compute the economical SVD with optimal complexity and store
  // the results in U, s and V.
  // START

  // Computation of \cor{economical QR-decompositions} of B-factors
  // incurs a computational effort of
  // $\cob{O(nk^2)}$
  Eigen::HouseholderQR<Eigen::MatrixXd> qrB1(B1);  // \Label{cmr:a}
  const Eigen::MatrixXd Q1{qrB1.householderQ() *
                           Eigen::MatrixXd::Identity(n, k)};  // \Label{cmr:1}
  Eigen::HouseholderQR<Eigen::MatrixXd> qrB2(B2);             // \Label{cmr:b}
  const Eigen::MatrixXd Q2{qrB2.householderQ() *
                           Eigen::MatrixXd::Identity(n, k)};  // \Label{cmr:2}
  // Computation of the matrix $\cob{\VZ}$ as defined in \prbeqref{eq:Zdef}
  // Uses \eigen's matrix concatenation operator

  Eigen::MatrixXd Z(m, 2 * k);
  Z << (A1 * (qrB1.matrixQR()
                  .block(0, 0, k, k)
                  .template triangularView<Eigen::Upper>())
                 .transpose()),
      (A2 * (qrB2.matrixQR()
                 .block(0, 0, k, k)
                 .template triangularView<Eigen::Upper>())
                .transpose());

  // Economical SVD of $\cob{\VZ}$
  Eigen::JacobiSVD<Eigen::MatrixXd> svdZ(
      Z, Eigen::ComputeThinU | Eigen::ComputeThinV);  // \Label{cmr:3}
  // Extract desired factors of final economical SVD of $\cob{\VX}$
  U = svdZ.matrixU();         // left factor of eco-SVD of $\cob{\VX}$
  s = svdZ.singularValues();  // singular values in a vector
  // Build right factor of final economical SVD
  Eigen::MatrixXd VZ{svdZ.matrixV()};
  V.block(0, 0, n, 2 * k) = Q1 * VZ.block(0, 0, k, 2 * k);  // \Label{cmr:4}
  V.block(n, 0, n, 2 * k) = Q2 * VZ.block(k, 0, k, 2 * k);  // \Label{cmr:5}
  // END

  return {U, s, V};
}
/* SAM_LISTING_END_1 */

/**
 * @brief Test for correct implementation of eco_svd_concat
 *
 * @param A1 matrix 1, left factor of low-rank factorization
 * @param B1 matrix 1, right factor of low-rank factorization
 * @param A2 matrix 2, left factor of low-rank factorization
 * @param B2 matrix 2, right factor of low-rank factorization
 * @return true if eco_svd_concat is correctly implemented
 * @return false otherwise
 */
/* SAM_LISTING_BEGIN_2 */
bool test_eco_svd_concat(const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
                         const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2) {
  auto [U, s, V] = eco_svd_concat(A1, B1, A2, B2);
  // Obtain matrix dimensions and perform consistency checks
  const unsigned int m = A1.rows();
  const unsigned int n = B1.rows();
  const unsigned int k = A1.cols();
  bool ok = ((U.rows() == m) && (U.cols() == 2 * k) && (V.rows() == 2 * n) &&
             (V.cols() == 2 * k) && (s.size() == 2 * k));
  // Define a suitable tolerance
  const double tau = 100 * m * n * std::numeric_limits<double>::epsilon();
  // Ensure that the singular values are non-negative and ordered
  ok = ok && (s[0] >= 0.0);
  for (unsigned int i = 1; i < 2 * k; ++i) {
    ok = ok && (s[i - 1] >= s[i]) && (s[i] >= 0.0);
  }

  // TODO: (3-15.c) Check orthonormality and correctness.
  // START
  // Check whether the columns of U are orthonormal
  ok = ok &&
       ((U.transpose() * U - Eigen::MatrixXd::Identity(2 * k, 2 * k)).norm() <=
        tau);
  // Check whether the columns of V are orthonormal
  ok = ok &&
       ((V.transpose() * V - Eigen::MatrixXd::Identity(2 * k, 2 * k)).norm() <=
        tau);
  // Check correctness of matrix factorization
  const Eigen::MatrixXd P{U * s.asDiagonal() * V.transpose()};
  const double ref = tau * s.norm();
  ok = ok && ((A1 * B1.transpose() - P.block(0, 0, m, n)).norm() <= ref);
  ok = ok && ((A2 * B2.transpose() - P.block(0, n, m, n)).norm() <= ref);
  // END

  return ok;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Best low-rank approximation of concatenated low-rank matrices
 *
 * @param A1 matrix 1, left factor of low-rank factorization
 * @param B1 matrix 1, right factor of low-rank factorization
 * @param A2 matrix 2, left factor of low-rank factorization
 * @param B2 matrix 2, right factor of low-rank factorization
 * @param tol tolerance
 * @return std::pair<Eigen::MatrixXd, Eigen::MatrixXd> low-rank factorization of
 * solution matrix
 */
/* SAM_LISTING_BEGIN_3 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> concat_low_rank_best(
    const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
    const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2, double tol) {
  assert((tol > 0) && (tol < 1.0));
  // Obtain \cor{economical SVD} of $\cob{\VX}$, \prbautoref{sp:2}
  auto [U, s, V] = eco_svd_concat(A1, B1, A2, B2);
  unsigned int L = s.size();
  Eigen::MatrixXd A_M = Eigen::MatrixXd::Zero(L, L);  // return dummy
  Eigen::MatrixXd B_M = Eigen::MatrixXd::Zero(L, L);  // return dummy

  // TODO: (3-15.d) Compute best low-rank approximation and return
  // factorization.
  // START

  // Find minimal rank $\cob{r}$, see
  // \prbeqref{eq:r}
  unsigned int r;
  for (r = 1; (r < L) && (s[r] / s[0] > tol); ++r)
    ;
  // Rescale columns of U
  for (unsigned int j = 0; j < r; ++j) U.col(j) *= s[j];
  // Extract sub-matrices
  A_M = U.leftCols(r);
  B_M = V.leftCols(r);
  // END
  return {A_M, B_M};
}
/* SAM_LISTING_END_3 */

#endif