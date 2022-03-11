/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#define _USE_MATH_DEFINES

// Function computing the factors of the economical SVD of a matrix created by
// concatenating two equal-size matrices given in low-rank factorization
/* SAM_LISTING_BEGIN_1 */
std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd>
eco_svd_concat(const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
               const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2) {
  // Matrix dimensions and consistency checks
  const int m = A1.rows();
  assert(A2.rows() == m);
  const int n = B1.rows();
  assert(B2.rows() == n);
  const int k = A1.cols();
  assert((k == A2.cols()) && (k == B1.cols()) && (k == B2.cols()));
  // Matrices/vector to be returned by the function
  Eigen::MatrixXd U(m, 2 * k);
  Eigen::VectorXd s(2 * k);
  Eigen::MatrixXd V(2 * n, 2 * k);
  // START: Student solution
  // Computation of \cor{economical QR-decompositions} of B-factors
  // Incurs a computational effort of $\cob{O(nk^2)}$
  Eigen::HouseholderQR<Eigen::MatrixXd> qrB1(B1); // \Label{cmr:a}
  const Eigen::MatrixXd Q1{qrB1.householderQ() *
                           Eigen::MatrixXd::Identity(n, k)}; // \Label{cmr:1}
  Eigen::HouseholderQR<Eigen::MatrixXd> qrB2(B2);            // \Label{cmr:b}
  const Eigen::MatrixXd Q2{qrB2.householderQ() *
                           Eigen::MatrixXd::Identity(n, k)}; // \Label{cmr:2}
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
      Z, Eigen::ComputeThinU | Eigen::ComputeThinV); // \Label{cmr:3}
  // Extract desired factors of final economical SVD of $\cob{\VX}$
  U = svdZ.matrixU();        // left factor of eco-SVD of $\cob{\VX}$
  s = svdZ.singularValues(); // singular values in a vector
  // Build right factor of final economical SVD
  Eigen::MatrixXd VZ{svdZ.matrixV()};
  V.block(0, 0, n, 2 * k) = Q1 * VZ.block(0, 0, k, 2 * k); // \Label{cmr:4}
  V.block(n, 0, n, 2 * k) = Q2 * VZ.block(k, 0, k, 2 * k); // \Label{cmr:5}
  // END: Student solution
  return {U, s, V};
}
/* SAM_LISTING_END_1 */

// Test for correct implementation of eco_svd_concat
/* SAM_LISTING_BEGIN_2 */
bool test_eco_svd_concat(const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
                        const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2) {
  auto [U, s, V] = eco_svd_concat(A1, B1, A2, B2);
  // Obtain matrix dimensions and perform consistency checks
  const int m = A1.rows();
  const int n = B1.rows();
  const int k = A1.cols();
  bool ok = ((U.rows() == m) && (U.cols() == 2 * k) && (V.rows() == 2 * n) &&
            (V.cols() == 2 * k) && (s.size() == 2 * k));
  // Define a suitable tolerance
  const double tau = 100 * m * n * std::numeric_limits<double>::epsilon();
  // Ensure that the singular values are non-negative and ordered
  ok = ok && (s[0] >= 0.0);
  for (int i = 1; i < 2 * k; ++i) {
    ok = ok && (s[i - 1] >= s[i]) && (s[i] >= 0.0);
  }
  // START: Student solution
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
  // END Student solution
  return ok;
}
/* SAM_LISTING_END_2 */

// Best low-rank approximation of concatenated low-rank matrices
/* SAM_LISTING_BEGIN_3 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
concat_low_rank_best(const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
                    const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2,
                    double tol) {
  assert((tol > 0) && (tol < 1.0));
  // Obtain \cor{economical SVD} of $\cob{\VX}$, \prbautoref{sp:2}
  auto [U, s, V] = eco_svd_concat(A1, B1, A2, B2);
  int L = s.size();
  // START student solution
  // Find minimal rank $\cob{r}$, see \prbeqref{eq:r}
  int r;
  for (r = 1; (r < L) && (s[r] / s[0] > tol); ++r)
    ;
  // Rescale columns of U
  for (int j = 0; j < r; ++j)
    U.col(j) *= s[j];
  // Extract sub-matrices
  return {U.leftCols(r), V.leftCols(r)};
  // END student solution
}
/* SAM_LISTING_END_3 */
