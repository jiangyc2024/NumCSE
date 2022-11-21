#ifndef EXPONENTIALINTEGRATOR_H_
#define EXPONENTIALINTEGRATOR_H_

/**
 * @file exponentialintegrator.cc
 * @brief NPDE homework ExponentialIntegrator code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

// Function $\phi$ used in the Exponential Euler
// single step method for an autonomous ODE.
Eigen::MatrixXd phim(const Eigen::MatrixXd &Z) {
  int n = Z.cols();
  assert(n == Z.rows() && "Matrix must be square.");
  Eigen::MatrixXd C(2 * n, 2 * n);
  C << Z, Eigen::MatrixXd::Identity(n, n), Eigen::MatrixXd::Zero(n, 2 * n);
  return C.exp().block(0, n, n, n);
}

// Calculate a single step of the exponential Euler method.
/* SAM_LISTING_BEGIN_0 */
template <class Function, class Jacobian>
Eigen::VectorXd exponentialEulerStep(const Eigen::VectorXd &y0, Function &&f,
                                     Jacobian &&df, double h) {
  // TO DO: 12-5.d
  // START
  return y0 + h * phim(h * df(y0)) * f(y0);
  // END
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void testExpEulerLogODE() {
  // TO DO: 12-5.e
  // START
  // Final time
  double T = 1.0;
  // Initial value
  Eigen::VectorXd y0(1);
  y0 << 0.1;
  // Function and Jacobian and exact solution
  auto f = [](const Eigen::VectorXd &y) { return y(0) * (1.0 - y(0)); };
  auto df = [](const Eigen::VectorXd &y) {
    Eigen::MatrixXd dfy(1, 1);
    dfy << 1.0 - 2.0 * y(0);
    return dfy;
  };
  double exactyT = y0(0) / (y0(0) + (1.0 - y0(0)) * std::exp(-T));

  // Container for errors
  std::vector<double> error(15);

  // Test many step sizes
  for (int j = 0; j < 15; ++j) {
    int M = std::pow(2, j + 1);
    Eigen::VectorXd y = y0;
    double h = T / M;
    for (int k = 0; k < M; ++k) {
      y = exponentialEulerStep(y, f, df, h);
    }

    error[j] = std::abs(y(0) - exactyT);
    std::cout << std::left << std::setfill(' ') << std::setw(3)
              << "M = " << std::setw(7) << M << std::setw(8)
              << "Error = " << std::setw(13) << error[j];
    if (j > 0) {
      std::cout << std::left << std::setfill(' ') << std::setw(10)
                << "Approximated order = " << std::log2(error[j - 1] / error[j])
                << std::endl;
    } else
      std::cout << std::endl;
  }
  // END
}
/* SAM_LISTING_END_1 */

#endif // #ifndef EXPONENTIALINTEGRATOR_H_
