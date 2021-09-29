#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "solvepermb.hpp"
#include "timer.h"

int main() {
  constexpr unsigned int n = 9;
  // Compute with both solvers
  std::cout << "--> Check that the solvers are correct" << std::endl;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n);
  Eigen::VectorXd b = Eigen::VectorXd::Random(n);
  Eigen::MatrixXd Xi, Xr, X;

  solvpermb(A, b, Xi);
  solvpermb_on3(A, b, Xr);
  std::cout << "Error = " << (Xi - Xr).norm() << std::endl;

  // Compute runtimes of different solvers
  std::cout << "--> Runtime comparison of naive solver vs reusing LU"
            << std::endl;
  constexpr unsigned int repeats = 3;

  // Header
  std::cout << std::setw(20) << "n" << std::setw(20) << "time no reuse [s]"
            << std::setw(20) << "time reuse [s]" << std::endl;

  // Loop over matrix size
  for (unsigned int p = 2; p <= 7; ++p) {
    // Timers
    Timer tm_naive, tm_reuseLU;
    unsigned int n = pow(2, p);

    // Repeat test many times
    for (unsigned int r = 0; r < repeats; ++r) {
      A = Eigen::MatrixXd::Random(n, n);
      b = Eigen::VectorXd::Random(n);

      // Compute runtime with naive solver
      tm_naive.start();
      solvpermb(A, b, X);
      tm_naive.stop();
      // Compute runtime reusing LU factorisation
      tm_reuseLU.start();
      solvpermb_on3(A, b, X);
      tm_reuseLU.stop();
    }
    // Print runtimes
    std::cout << std::setw(20) << n << std::scientific << std::setprecision(3)
              << std::setw(20) << tm_naive.min() << std::setw(20)
              << tm_reuseLU.min() << std::endl;
  }
}