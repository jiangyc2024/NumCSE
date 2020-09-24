///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include "spy.hpp"

int main() {
  /* SAM_LISTING_BEGIN_0 */
  // Study of fill-in with LU-factorization due to pivoting
  MatrixXd A(11, 11);
  A.setZero();
  A.diagonal() = VectorXd::LinSpaced(11, 1, 11).cwiseInverse();
  A.col(10).setConstant(2);
  A.row(10).setConstant(2);
  auto solver = A.lu();
  MatrixXd L = MatrixXd::Identity(11, 11);
  L += solver.matrixLU().triangularView<StrictlyLower>();
  MatrixXd U = solver.matrixLU().triangularView<Upper>();
  // Plotting
  spy(A, "Arrow matrix A", "fillinpivotA.eps");
  spy(L, "L factor", "fillinpivotL.eps");
  spy(U, "U factor", "fillinpivotU.eps");
  std::cout << A << std::endl;
  /* SAM_LISTING_END_0 */
  return 0;
}
