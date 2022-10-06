#ifndef COMPACTSTORAGEQR_HPP
#define COMPACTSTORAGEQR_HPP

/*
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

/**
 * @brief wrapper for compactly stored QR-factorization
 *
 */
/* SAM_LISTING_BEGIN_1 */
class CompactStorageQR {
 public:
  /**
   * @brief Construct a new CompactStorageQR object
   *
   * @param data raw array which has to be persistent
   * @param n number of rows/columns of encoded matrix
   */
  CompactStorageQR(const double *data, unsigned int n) : n_(n), M_(data, n, n) {
    for (unsigned int k = 0; k < n_ - 1; ++k) {
      double sn{0};
      for (unsigned int j = k + 1; j < n_; ++j) {
        const double t{M_(j, k)};
        sn += t * t;
      }
      if (sn > 1.0) {
        throw std::runtime_error(
            "CompactStorageQR: Illegal subdiagonal column norm!");
      }
    }
  }

  /**
   * @brief Determinant of the matrix
   *
   * @return double the determinant
   */
  double det() const;

  /**
   * @brief  Right multiplication with another matrix
   *
   * @param X matrix to multiply
   * @return Eigen::MatrixXd product
   */
  Eigen::MatrixXd matmult(const Eigen::MatrixXd &X) const;

  /**
   * @brief Solution of linear systems of equations
   *
   * @param B right hand side
   * @return Eigen::MatrixXd solution matrix
   */
  Eigen::MatrixXd solve(const Eigen::MatrixXd &B) const;

 private:
  unsigned int n_;                       // Matrix dimensions
  Eigen::Map<const Eigen::MatrixXd> M_;  // Raw data wrapped into a matrix
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::MatrixXd CompactStorageQR::matmult(const Eigen::MatrixXd &X) const {
  assert((X.rows() == n_) && "Wrong size of matrix factor");
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(n_, X.rows());

  // TODO: (3-16.a) Implement matrix multiplication with computational cost at
  // most O(n^2k)
  // START

  // END

  return Y;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd CompactStorageQR::solve(const Eigen::MatrixXd &B) const {
  assert((B.rows() == n_) && "Wrong size of right-hand side");
  Eigen::MatrixXd X{B};

  // TODO: (3-16.b) Solve the linear system of equations given by AX = B with
  // complexity O(n^2k).
  // START

  // END

  return X;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
double CompactStorageQR::det() const {
  double d = 0.;

  // TODO: (3-16.c) Compute the determinant of A.
  // START

  // END

  return d;
}
/* SAM_LISTING_END_4 */

#endif