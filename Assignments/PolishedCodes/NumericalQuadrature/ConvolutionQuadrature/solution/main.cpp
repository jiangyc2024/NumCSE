/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <unsupported/Eigen/FFT>
#include <vector>

#include "convolutionquadrature.hpp"

int main() {
  std::cout << "NumCSE code for computation of convolution quadrature weights"
            << std::endl;

  // Test case: Laplace transform is a square root
  auto F = [](std::complex<double> s) -> std::complex<double> {
    return std::sqrt(s);
  };
  {
    std::cout << "Convolution quadrature weights for square-root function"
              << std::endl;
    std::cout << "\n RESULT weights w = "
              << compute_cq_weights(F, 20, 0.1).transpose().real() << std::endl;
    std::cout << "REFERENCE weights = "
              << compute_cq_weights_sqrt(20, 0.1).transpose() << std::endl;
  }
  {
    std::cout << "N-convergence of convolution quadrature weights" << std::endl;
    constexpr double r = 0.99;
    constexpr double tau = 0.1;
    for (unsigned int N = 3; N < 25; ++N) {
      Eigen::VectorXd w_comp = compute_cq_weights(F, N, tau, r).real();
      Eigen::VectorXd w_ref = compute_cq_weights_sqrt(N, tau);
      std::cout << std::setw(15) << std::scientific
                << std::abs(w_comp[1] - w_ref[1]) << " & " << std::setw(15)
                << std::scientific << std::abs(w_comp[2] - w_ref[2]) << " & "
                << std::setw(15) << std::scientific
                << std::abs(w_comp[3] - w_ref[3]) << " \\\\" << std::endl;
    }
  }
  return 0;
}
