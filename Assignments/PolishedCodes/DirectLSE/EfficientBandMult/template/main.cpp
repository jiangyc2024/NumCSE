#include <Eigen/Dense>
#include <iostream>

#include "efficientbandmult.hpp"

int main() {
  constexpr unsigned int n = 9;
  // Compute with all three solvers
  std::cout << "*** Check that the solvers are correct" << std::endl;
  Eigen::VectorXd a = Eigen::VectorXd::Random(n - 1);
  Eigen::VectorXd b =
      Eigen::VectorXd::Zero(n - 2);  // All 0s for upper diagonal structure
  Eigen::VectorXd y = Eigen::VectorXd::Random(n);
  Eigen::VectorXd x;

  std::cout << "Original: " << y << std::endl;

  // Compute $y = A*A^{-1}*y$
  solvelseAupper(a, y, x);
  multAx(a, b, x, y);

  // Should be the same as before
  std::cout << "Upper: " << y << std::endl;

  // Random $b$ for own and Eigen-based solver
  b = Eigen::VectorXd::Random(n - 2);

  // Compute $y = A*A^{-1}*y$
  solvelseA(a, b, y, x);
  multAx(a, b, x, y);

  // Should be the same as before
  std::cout << "Own: " << y << std::endl;

  // Compute $y = A*A^{-1}*y$
  solvelseAEigen(a, b, y, x);
  multAx(a, b, x, y);

  // Should be the same as before
  std::cout << "Eigen: " << y << std::endl;
}