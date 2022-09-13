/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#include <Eigen/Dense>
#include <Eigen/QR>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

// Wrapper class for compactly stored QR-factorization
class CompactStorageQR {
 public:
  // Constructor taking raw data array, which has to be persistent
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
  // Determinant of the matrix stored in this object
  double det() const;
  // Right multiplication with another matrix
  Eigen::MatrixXd matmult(const Eigen::MatrixXd &X) const;
  // Solution of linear systems of equations
  Eigen::MatrixXd solve(const Eigen::MatrixXd &B) const;
 private:
  unsigned int n_;                       // Matrix dimensions
  Eigen::Map<const Eigen::MatrixXd> M_;  // Raw data wrapped into a matrix
};

Eigen::MatrixXd CompactStorageQR::matmult(const Eigen::MatrixXd &X) const {
  assert((X.rows() == n_) && "Wrong size of matrix factor");
  Eigen::MatrixXd Y(n_, X.rows());
  // START Student code
  // Multiplication with R-factor
  Y = M_.template triangularView<Eigen::Upper>() * X;
  // Multiplication with Q-factors
  Eigen::VectorXd u{Eigen::VectorXd::Zero(n_)};
  for (int k = n_ - 2; k >= 0; --k) {
    u.tail(n_ - k - 1) = M_.block(k + 1, k, n_ - k - 1, 1);
    u[k] = std::sqrt(1.0 - u.squaredNorm());
    // std::cout<<"output u"<< u << std::endl;
    Y -= 2 * u * (u.transpose() * Y);
  }
  // END Student code
  // std::cout<<"output Y"<< Y << std::endl;
  return Y;
}

Eigen::MatrixXd CompactStorageQR::solve(const Eigen::MatrixXd &B) const {
  assert((B.rows() == n_) && "Wrong size of right-hand side");
  Eigen::MatrixXd X{B};
  // START student code
  // Check invertibility by inspecting diagonal elements of R-factor
  const double atol = 1.0E-16;
  const double rtol = 1.0E-8;
  double maxrii = std::abs(M_(0, 0));
  for (unsigned int k = 1; k < n_; ++k) {
    const double tmp = std::abs(M_(k, k));
    if (maxrii < tmp) {
      maxrii = tmp;
    }
  }
  for (unsigned int k = 1; k < n_; ++k) {
    const double rii = std::abs(M_(k, k));
    if ((rii < rtol * maxrii) || (rii < atol)) {
      throw std::runtime_error("CompactStorageQR::solve: matrix not regular");
    }
  }
  // Multiplication with transposed Q-factors
  Eigen::VectorXd u(n_);
  for (unsigned int k = 0; k < n_ - 1; ++k) {
    double sn{0};
    for (unsigned int j = k + 1; j < n_; ++j) {
      u[j] = M_(j, k);
      sn += u[j] * u[j];
    }
    u[k] = std::sqrt(1.0 - sn);
    std::cout<<"output here"<< u << std::endl;
    // Multiplication with current Q-factor
    X -= 2 * u * (u.transpose() * X);
    u[k] = 0.0;
  }
  // Solve upper triangular linear system $\cob{\VR\VX=\VQ^{\top}\VB}$
  X = M_.template triangularView<Eigen::Upper>().solve(X);
  // END student code
  return X;
}

double CompactStorageQR::det() const {
  double d;
  // START Student code
  // The determinant is just the product of the diagonal entries of the
  // R-factor, because the determinant of the Q-factor is 1.
  d = 1.0;
  for (unsigned int j = 0; j < n_; ++j) {
    d *= M_(j, j);
  }
  // Sign correction, because Househiolder transformations
  // have determinant -1 and n-1 of them are inviolved.
  d *= ((n_%2 == 0)?-1.0:1.0);
  // END Student code
  return d;
}
