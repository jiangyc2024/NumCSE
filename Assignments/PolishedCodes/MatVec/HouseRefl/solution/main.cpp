#include <Eigen/Dense>
#include <ctime>
#include <iostream>

#include "houserefl.hpp"

int main() {
  // Parameters
  constexpr int n = 6;
  constexpr double eps = 1e-10;

  // Allocate memory
  Eigen::VectorXd v(n);
  Eigen::MatrixXd Z(n, n - 1);

  // Initialize random number generator
  srand((unsigned int)time(0));

  // Fill v with random values
  v = Eigen::VectorXd::Random(n);

  // Get ONB of the space orthogonal to v
  std::cout << "Running algorithm for v=\n" << v << "\n\n";
  houserefl(v, Z);
  std::cout << "Result is Z=\n" << Z << "\n\n";

  // Compute error
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(n - 1, n - 1);
  double error = (Z.transpose() * Z - Id).norm();

  // Check correctness
  std::cout << "Error = " << error << "\n\n";

  if (error < eps) {
    std::cout << "Your implementation seems to be correct.\n";
  } else {
    std::cout << "The error is bigger than 10e-10.\n";
    std::cout << "Your implementation seems to be wrong.\n";
  }
}
