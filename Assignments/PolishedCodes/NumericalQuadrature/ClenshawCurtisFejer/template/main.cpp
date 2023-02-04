#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <unsupported/Eigen/FFT>

#include "clenshawcurtisfejer.hpp"

int main() {
  std::cout << "C++ code for NumCSE exam problem on CCF quadrature rules"
            << std::endl;

  // Printing nodes and weights of quadrature rules
  std::cout << "List quadrature weights and nodes" << std::endl;
  for (unsigned int n = 1; n < 6; n++) {
    const CCFQuadRule ccfqr(n);
    const CCFQuadRule_Fast ccfqr_fast(n);
    std::cout << "n = " << n << " : " << std::endl
              << "\tNodes = " << ccfqr.nodes().transpose() << std::endl
              << "\tWeights = " << ccfqr.weights().transpose() << std::endl
              << "\tFast weights = " << ccfqr_fast.weights().transpose()
              << std::endl;
  }

  return 0;
}
