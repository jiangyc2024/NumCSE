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
  return y0;
  // END
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void testExpEulerLogODE() {
  // TO DO: 12-5.e
  // START

  // END
}
/* SAM_LISTING_END_1 */

#endif // #ifndef EXPONENTIALINTEGRATOR_H_
