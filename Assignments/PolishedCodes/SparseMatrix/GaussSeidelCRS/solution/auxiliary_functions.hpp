#ifndef AUXFUNCSHPP
#define AUXFUNCSHPP
/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "gaussseidelcrs.hpp"

// One step of a Gauss-Seidel iteration for the linear system of equations Ax=b
// with A given as some Eigen matrix type

/* SAM_LISTING_BEGIN_1 */
template <typename MATRIX>
bool GaussSeidelstep_generic(const MATRIX &A, const Eigen::VectorXd &b,
                             Eigen::VectorXd &x) {
  const unsigned int n = A.cols();
  assert(n == A.rows() && "Matrix must be square");
  assert(n == b.size() && "Vector length mismatch");
  assert(n == x.size() && "Vector b length mismatch");

  for (unsigned int i = 0; i < n; ++i) {
    if (A(i, i) == 0.0) return false;
    double s = b[i];
    for (unsigned int j = 0; j < i; ++j) s -= A(i, j) * x[j];
    for (unsigned int j = i + 1; j < n; ++j) s -= A(i, j) * x[j];
    x[i] = s / A(i, i);
  }
  return true;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_4 */
template <typename RECORDER = std::function<void(double, double)>>
bool GaussSeidel_iteration_II(
    const CRSMatrix &A, const Eigen::VectorXd &b, Eigen::VectorXd &x,
    double atol = 1.0E-8, double rtol = 1.0E-6, unsigned int maxit = 100,
    RECORDER &&rec = [](double, double) -> void {}) {
  assert(A.n == A.m && "Matrix must be square");
  assert(A.n == b.size() && "Vector b length mismatch");
  assert(A.n == x.size() && "Vector x length mismatch");

  // Main loop for Gauss-Seidel iteration
  Eigen::VectorXd x_old{x};
  for (unsigned int k = 0; k < maxit; ++k) {
    // TODO: (2-18.d; optional) Second function to implement the Gauss-Seidel
    // iteration, not required to implement.
    // START
    if (!GaussSeidelstep_crs(A, b, x)) return false;
    const double delta_x_norm = (x - x_old).norm();
    const double x_norm = x.norm();
    rec(x_norm, delta_x_norm);  // \Label{xgs:r}
    if ((delta_x_norm < atol) ||
        (delta_x_norm < rtol * x_norm)) {  // \Label{xgs:1}
      return true;
    }
    x_old = x;
    // END
  }
  // No termination
  return false;
}
/* SAM_LISTING_END_4 */

#endif
