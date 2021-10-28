#include <Eigen/Dense>
#include <iostream>

#include "blockLU.hpp"

int main() {
  constexpr unsigned int n = 10;
  Eigen::VectorXd v, u, b;
  u = v = Eigen::VectorXd::Random(n);
  b = Eigen::VectorXd::Random(n + 1);
  // upper triangular matrix
  Eigen::MatrixXd R(n, n);
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = i; j < n; ++j) {
      R(i, j) = rand();  // Bad RNG, but sufficient here
    }
  }
  R /= RAND_MAX;  // "norm" R for numerical stability
                  // Build matrix A for Eigensolver
  Eigen::MatrixXd A(n + 1, n + 1);
  A << R, v, u.transpose(), 0;

  const double error =
      (solve_LSE(R, v, u, b) - A.colPivHouseholderQr().solve(b)).norm();
  if (error > 1e-8) {
    std::cout << "solve_LSE() returns a different solution than Eigen"
              << std::endl;
  } else {
    std::cout << "solve_LSE() and Eigen get the same result" << std::endl;
  }
}