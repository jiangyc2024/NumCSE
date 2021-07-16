#include <Eigen/Dense>
#include <iostream>
#include <limits>

#include "choleskyQR.hpp"

int main() {
  constexpr size_t m = 3;
  constexpr size_t n = 2;

  Eigen::MatrixXd A(m, n);
  constexpr double epsilon = std::numeric_limits<double>::denorm_min();
  A << 3, 5, 1, 9, 7, 1;
  // A << 1, 1, 0.5 * epsilon, 0, 0, 0.5 * epsilon;
  std::cout << "A = " << std::endl << A << std::endl;

  Eigen::MatrixXd R, Q;
  CholeskyQR(A, R, Q);

  std::cout << "From Cholesky: R = " << std::endl << R << std::endl;
  std::cout << "From Cholesky: Q = " << std::endl << Q << std::endl;

  DirectQR(A, R, Q);

  std::cout << "Direct QR: R = " << std::endl << R << std::endl;
  std::cout << "Direct QR: Q = " << std::endl << Q << std::endl;
}