///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cassert>

#include <Eigen/Dense>

namespace lumult {


using Eigen::MatrixXd;

inline
/* SAM_LISTING_BEGIN_0 */
//! Multiplication of normalized lower/upper triangular matrices
MatrixXd lumult(const MatrixXd &L, const MatrixXd &U) {
  const Eigen::Index n = L.rows();
  assert(n == L.cols() && n == U.cols() && n == U.rows());
  MatrixXd A{MatrixXd::Zero(n, n)};
  for (Eigen::Index k = 0; k < n; ++k) {
    for (Eigen::Index j = k; j < n; ++j) {
      A(k, j) = U(k, j) + (L.block(k, 0, 1, k) * U.block(0, j, k, 1))(0, 0);
    }
    for (Eigen::Index i = k + 1; i < n; ++i) {
      A(i, k) =
          (L.block(i, 0, 1, k) * U.block(0, k, k, 1))(0, 0) + L(i, k) * U(k, k);
    }
  }
  return A;
}
/* SAM_LISTING_END_0 */


} //namespace lumult