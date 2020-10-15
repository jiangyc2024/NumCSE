#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>
#include <vector>

#include "totriplets.hpp"

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
SparseMatrix<double> initA(unsigned int n) {
  unsigned int rows = n * (n - 1) / 2;
  unsigned int cols = n - 1;
  SparseMatrix<double> A(rows, cols);
  // TODO: (4-10.b) Initialize the sparse coefficient matrix for
  // the distance fitting problem
  // START
  
  // END
  A.makeCompressed();
  return A;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
VectorXd solveExtendedNormalEquations(const MatrixXd &D) {
  VectorXd x;
  unsigned int n = D.cols();
  unsigned int m = n * (n - 1) / 2;
  // TODO: (4-10.c) Solve the extended normal equations with a sparse solver
  // that Eigen provides
  // START
  
  // END
  return x;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
VectorXd solveNormalEquations(const MatrixXd &D) {
  VectorXd x;
  // TODO: (4-10.e) Solve the normal equations by exploiting the
  // Sherman-Morrison-Woodbury formula, i.e. by using your results from the
  // previous subproblem
  // START
  
  // END
  return x;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_END_3 */
bool testNormalEquations(const MatrixXd &D) {
  // TODO: (4-10.f) Call your implementations of solveExtendedNormalEquations()
  // and solveNormalEquations() and return true, if their results agree.
  // START
  
  // END
  return false;
}
/* SAM_LISTING_END_3 */
