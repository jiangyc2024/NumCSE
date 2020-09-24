///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "matplotlibcpp.h"
#include "spdiags.hpp"
#include "spy.hpp"

using namespace std;
using namespace Eigen;

namespace plt = matplotlibcpp;

int main() {
  /* SAM_LISTING_BEGIN_0 */
  // Build matrix
  int n = 100;
  RowVectorXd diag_el(5);
  diag_el << -1, -1, 3, -1, -1;
  VectorXi diag_no(5);
  diag_no << -n, -1, 0, 1, n;
  MatrixXd B = diag_el.replicate(2 * n, 1);
  B(n - 1, 1) = 0;
  B(n, 3) = 0;  // delete elements
  // A custom function from the Utils folder
  SparseMatrix<double> A = spdiags(B, diag_no, 2 * n, 2 * n);
  // It is not possible to access the LU-factors in the case of
  // \eigen's LU-decomposition for sparse matrices.
  // Therefore we have to resort to the dense version.
  auto solver = MatrixXd(A).lu();
  MatrixXd L = MatrixXd::Identity(2 * n, 2 * n);
  L += solver.matrixLU().triangularView<StrictlyLower>();
  MatrixXd U = solver.matrixLU().triangularView<Upper>();
  // Plotting
  spy(A, "Sparse matrix", "sparseA_cpp.eps");
  spy(L, "Sparse matrix: L factor", "sparseL_cpp.eps");
  spy(U, "Sparse matrix: U factor", "sparseU_cpp.eps");
  /* SAM_LISTING_END_0 */

  return 0;
}
