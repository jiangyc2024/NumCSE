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
  // TO DO: Implement efficient matrix multiplication A*X where A is the matrix
  //        represented by the object CompactStorageQR
  // START 
  
  // END 
  return Y;
}

Eigen::MatrixXd CompactStorageQR::solve(const Eigen::MatrixXd &B) const {
  assert((B.rows() == n_) && "Wrong size of right-hand side");
  Eigen::MatrixXd X{B};
  // TO DO: Implement efficient function to solve A^-1*B where A is the matrix
  //        represented by the object CompactStorageQR
  // START 

  // END 
  return X;
}

double CompactStorageQR::det() const {
  double d;
  // TO DO: Implement a function to compute the determinant of the matrix
  //        represented by the object CompactStorageQR
  // START 

  // END 
  return d;
}
