#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#include "kron.hpp"
#include "timer.h"

int main() {
  // Testing correctness of Kron
  Eigen::MatrixXd A(2, 2);
  A << 1, 2, 3, 4;
  Eigen::MatrixXd B(2, 2);
  B << 5, 6, 7, 8;
  Eigen::MatrixXd C;
  // Fix the random seed for reproducibility
  srand(5);
  Eigen::VectorXd x = Eigen::VectorXd::Random(4);
  Eigen::VectorXd y;

  kron(A, B, C);
  std::cout << std::setprecision(4) << "kron(A,B) = " << std::endl
            << C << std::endl;

  std::cout << "\nEnter \"0\" to test all functions.\n"
            << "Enter \"1\" to only test kron().\n"
            << "Enter \"2\" to only test kron_mult().\n"
            << "Enter \"3\" to only test kron_reshape().\n"
            << "Enter \"4\" to only test kron_runtime().\n";

  int ans = 0;
  std::cin >> ans;
  switch (ans) {
    case 0:
      y = C * x;
      std::cout << "Using kron: y =\n" << y << std::endl;
      kron_mult(A, B, x, y);
      std::cout << "Using kron_mult: y =\n" << y << std::endl;
      kron_reshape(A, B, x, y);
      std::cout << "Using kron_reshape: y =\n" << y << std::endl;
      kron_runtime();
      break;
    case 1:
      y = C * x;
      std::cout << "Using kron: y =\n" << y << std::endl;
      break;
    case 2:
      kron_mult(A, B, x, y);
      std::cout << "Using kron_mult: y =\n" << y << std::endl;
      break;
    case 3:
      kron_reshape(A, B, x, y);
      std::cout << "Using kron_reshape: y =\n" << y << std::endl;
      break;
    case 4:
      kron_runtime();
      break;
  }
}
