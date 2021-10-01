#include <Eigen/Dense>
#include <iostream>

#include "blocklsepiv.hpp"

int main() {
  constexpr unsigned int n = 5;
  Eigen::VectorXd d1, d2, c, b(2 * n);
  d1 = Eigen::VectorXd::LinSpaced(n, 1, n);
  d2 = -d1;
  c = Eigen::VectorXd::Ones(n);
  b << d1, d1;

  std::cout << "--> multA" << std::endl;
  Eigen::VectorXd y = multA(d1, d2, c, b);
  std::cout << y << std::endl;

  std::cout << "--> solveA" << std::endl;
  y = solveA(d1, d2, c, b);
  std::cout << y << std::endl;

  std::cout << "--> numerical experiment" << std::endl;
  numericalExperiment();
}