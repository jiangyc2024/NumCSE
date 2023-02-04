#ifndef GAUSSSEIDELCRSHPP
#define GAUSSSEIDELCRSHPP
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

  // TODO: (2-18.b) Implement a single step of the Gauss-Seidel iteration with
  // the system matrix in CRS format.
  // START

  // END
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
    // TODO: (2-18.d) Implement the Gauss-Seidel iteration with a
    // correction-based termination criterion.
    // START

    // END
  }
  // No termination
  return false;
}
/* SAM_LISTING_END_5 */

#endif
