/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#ifdef NICEBACKTRACE
#include "backtrace.hpp"
#endif

#include <Eigen/Dense>
#include <Eigen/QR>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include "compactstorageqr.hpp"

#ifdef NICEBACKTRACE
#include "backtrace.hpp"
#endif

int main(int /*argc*/, char ** /*argv*/) {
  const int n = 10;
  // Generate test matrix
  Eigen::MatrixXd A(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      A(i, j) = 1.0/(2*i+j+1.0);
    }
  }
  std::cout << "Matrix A used to construct the CompactStorageQR object A_qr. A = \n \n" << A << std::endl;

  // Using member functions from the class CompactStorageQR
  CompactStorageQR A_qr(A.data(), n);
  
  // matmult
  Eigen::MatrixXd B = A_qr.matmult(Eigen::MatrixXd::Identity(n, n));
  std::cout << "\nMatrix B = A_qr * Id = \n \n" << B << std::endl;
  
  // solve
  Eigen::MatrixXd B_inv = A_qr.solve(Eigen::MatrixXd::Identity(n,n));
  std::cout << "\nMatrix B^-1 = \n \n" << B_inv << std::endl;
  
  // determinant
  std::cout << "\nDeterminant(A_qr) = " << A_qr.det() << std::endl;

  return 0;
}
