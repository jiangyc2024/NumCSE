#include <Eigen/Dense>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>

#include "matPow.hpp"

int main() {
  // Check/Test with provided, complex, matrix
  constexpr unsigned int n = 3;  // size of matrix
  constexpr unsigned int k = 9;  // power

  Eigen::MatrixXcd A = construct_matrix(n);

  // Output results
  Eigen::MatrixXcd Eigen_A_pow_k = A.pow(k);
  std::cout << "A = " << std::endl
            << A << std::endl
            << "-> Eigen:" << std::endl
            << "A^" << k << " = " << std::endl
            << Eigen_A_pow_k << std::endl;
  Eigen::MatrixXcd A_pow_k = matPow(A, k);
  std::cout << "-> Ours:" << std::endl
            << "A^" << k << " = " << std::endl
            << A_pow_k << std::endl;

  // Use this to test your implementation, error should be very small (machine
  // precision)
  std::cout << "-- Error of our implementation VS Eigen: "
            << (A_pow_k - Eigen_A_pow_k).norm() << std::endl;
  tabulateRuntime(n);
}