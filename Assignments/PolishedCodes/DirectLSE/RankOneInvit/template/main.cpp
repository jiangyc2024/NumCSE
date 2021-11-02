#include <Eigen/Dense>
#include <iostream>

#include "rankoneinvit.hpp"

int main() {
  // Initialization
  std::srand(21);
  constexpr double tol = 1e-3;
  double lmin, lmin_fast;
  constexpr unsigned int n = 10;

  // Compute with both implementations
  std::cout << "\n*** Test of both implementations\n";
  Eigen::VectorXd d = Eigen::VectorXd::Random(n);
  lmin = rankoneinvit(d, tol);
  lmin_fast = rankoneinvit_fast(d, tol);
  std::cout << "Slow implementation: " << lmin << std::endl;
  std::cout << "Fast implementation: " << lmin_fast << std::endl;
  std::cout << "Difference: " << (lmin_fast - lmin) << std::endl;

  std::cout << "\nEnter \"1\" to test rankoneinvit_runtime().\n";
  int ans = 0;
  std::cin >> ans;
  if (ans == 1) {
    // Compare runtimes of different implementations of rankoneinvit
    std::cout << "*** Runtime comparison of the two implementations"
              << std::endl;
    rankoneinvit_runtime();
  }
}