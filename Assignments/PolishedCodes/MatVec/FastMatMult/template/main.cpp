#include <Eigen/Dense>
#include <iostream>

#include "strassen.hpp"

int main() {
  constexpr double tolerance = 1e-9;
  constexpr int n = 128;
  std::srand(5);
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n);

  std::cout
      << "Works for identity matrix: "
      << ((strassenMatMult(A, Eigen::MatrixXd::Identity(n, n)) - A).norm() <
          tolerance)
      << " "
      << ((strassenMatMult(Eigen::MatrixXd::Identity(n, n), A) - A).norm() <
          tolerance)
      << std::endl;

  std::cout << "Works for two random matrices: "
            << (test_strassen() < tolerance) << std::endl;

  std::cout << "\nOutput of test_strassen()" << test_strassen() << std::endl;

  time_strassen();
}
