/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */
#ifndef GAUSSSEIDELCRSHPP
#define GAUSSSEIDELCRSHPP

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

#define _USE_MATH_DEFINES

struct CRSMatrix {
  unsigned int m;                     // number of rows
  unsigned int n;                     // number of columns
  std::vector<double> val;            // value array
  std::vector<unsigned int> col_ind;  // same length as value array
  std::vector<unsigned int> row_ptr;  // length m+1, row_ptr[m] == val.size()
};

// One step of a Gauss-Seidel iteration for the linear system of equations Ax=b
// with A given in a raw CRS format
/* SAM_LISTING_BEGIN_2 */
bool GaussSeidelstep_crs(const CRSMatrix &A, const Eigen::VectorXd &b,
                         Eigen::VectorXd &x) {
  assert(A.n == A.m && "Matrix must be square");
  assert(A.n == b.size() && "Vector b length mismatch");
  assert(A.n == x.size() && "Vector x length mismatch");

  // START student code
  // Outer loop over rows of the matrix
  for (unsigned int i = 0; i < A.n; ++i) {
    double Aii = 0.0;
    double s = b[i];
    // Inner summation loop over non-zero entries of i-th row.
    // Skip diagonal entry in the summation and store is separately.
    for (unsigned int l = A.row_ptr[i]; l < A.row_ptr[i + 1]; ++l) {
      const unsigned int j = A.col_ind[l];
      if (j != i) {
        s -= A.val[l] * x[j];
      } else {
  // Fetch diagonal entry of A.
        Aii = A.val[l];
      }
    }
    if (Aii != 0.0)
      x[i] = s / Aii;
    else
      return false;
  }
  // END student code
  return true;
}
/* SAM_LISTING_END_2 */

// Gauss-Seidel iteration for matrix A and vector b with correction-based
// termination
/* SAM_LISTING_BEGIN_5 */
bool GaussSeidel_iteration(const CRSMatrix &A, const Eigen::VectorXd &b,
                           Eigen::VectorXd &x, double atol = 1.0E-8,
                           double rtol = 1.0E-6, unsigned int maxit = 100) {
  assert(A.n == A.m && "Matrix must be square");
  assert(A.n == b.size() && "Vector b length mismatch");
  assert(A.n == x.size() && "Vector x length mismatch");

  // Main loop for Gauss-Seidel iteration
  for (unsigned int k = 0; k < maxit; ++k) {
    // START student code
    double x_normsq = 0.0;
    double delta_x_normsq = 0.0;
    for (unsigned int i = 0; i < A.n; ++i) {
      double Aii = 0.0;
      double s = b[i];
      for (unsigned int l = A.row_ptr[i]; l < A.row_ptr[i + 1]; ++l) {
        const unsigned int j = A.col_ind[l];
        if (j != i) {
          s -= A.val[l] * x[j];
        } else {
    // Fetch diagonal entry of A
          Aii = A.val[l];
        }
      }
      if (Aii != 0.0) {
        // Update iterate and compute Eucidean norm of current iterate
        // and of correction
        const double tmp = x[i];
        x[i] = s / Aii;
        delta_x_normsq += (x[i] - tmp) * (x[i] - tmp);
        x_normsq += x[i] * x[i];
      } else
        return false;
    }
    // Correction based termination cirterion
    if ((delta_x_normsq < atol * atol) ||
        (delta_x_normsq < rtol * rtol * x_normsq)) {
      return true;
    }
    // END student code
  }
  // No termination
  return false;
}
/* SAM_LISTING_END_5 */

#endif
