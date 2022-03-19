/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */
#ifndef AUXFUNCSHPP
#define AUXFUNCSHPP

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
template <typename MATRIX>
bool GaussSeidelstep_generic(const MATRIX &A, const Eigen::VectorXd &b,
                             Eigen::VectorXd &x) {
  const int n = A.cols();
  assert(n == A.rows() && "Matrix must be square");
  assert(n == b.size() && "Vector length mismatch");
  assert(n == x.size() && "Vector b length mismatch");

  for (int i = 0; i < n; ++i) {
    if (A(i, i) == 0.0) return false;
    double s = b[i];
    for (int j = 0; j < i; ++j) s -= A(i, j) * x[j];
    for (int j = i + 1; j < n; ++j) s -= A(i, j) * x[j];
    x[i] = s / A(i, i);
  }
  return true;
}
#endif
