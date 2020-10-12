#pragma once

#include "totriplets.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
SparseMatrix<double> initA(unsigned int n) {
  unsigned int rows = n * (n - 1) / 2;
  unsigned int cols = n - 1;
  SparseMatrix<double> A(rows, cols);
  // TO DO (0-1.b)
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
  // TO DO (0-1.c)
  // START
  
  // END
  return x;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
VectorXd solveNormalEquations(const MatrixXd &D) {
  VectorXd x;
  // TO DO (0-1.e)
  // START
  
  // END
  return x;
}

/* SAM_LISTING_END_2 */
