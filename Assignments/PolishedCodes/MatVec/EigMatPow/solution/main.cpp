#include <Eigen/Dense>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>

#include "getit.hpp"

int main() {
  // Some arbitrary data to test getit
  Eigen::MatrixXd A(4, 4);
  A << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
  Eigen::VectorXd x(4);
  x << 4, 5, 6, 7;
  constexpr unsigned int k = 9;

  // Testing the implementation with some matrix
  Eigen::VectorXd yg = getit(A, x, k);
  std::cout << "getit(A,x, k) = " << std::endl << yg << std::endl;

  // Checking that getit works
  Eigen::VectorXd yp = A.pow(k) * x;
  std::cout << "A^k x = " << std::endl << yp << std::endl;
  const double err = (yg - yp).norm() / yp.norm();
  std::cout << "Relative error = " << err << std::endl;
}