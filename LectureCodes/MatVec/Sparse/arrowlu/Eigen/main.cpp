///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "spy.hpp"
#include <Eigen/Dense>


using Eigen::MatrixXd;
using Eigen::StrictlyLower;
using Eigen::Upper;

//NOLINTBEGIN(bugprone-exception-escape)
int main() {
  /* SAM_LISTING_BEGIN_0 */
  // Build matrix
  MatrixXd A(11, 11);
  A.setIdentity();
  A.col(10).setOnes();
  A.row(10).setOnes();
  // A.reverseInPlace(); // used in\cref{ex:arrowmatrixlu}
  auto solver = A.lu();
  MatrixXd L = MatrixXd::Identity(11, 11);
  L += solver.matrixLU().triangularView<StrictlyLower>();
  const MatrixXd U = solver.matrixLU().triangularView<Upper>();
  const MatrixXd Ainv = A.inverse();
  // Plotting
  spy(A, "Pattern of A", "Apat_cpp.eps");
  spy(L, "Pattern of L", "Lpat_cpp.eps");
  spy(U, "Pattern of U", "Upat_cpp.eps");
  spy(Ainv, "Pattern of A^{-1}", "Ainvpat_cpp.eps");
  /* SAM_LISTING_END_0 */
  return 0;
}
//NOLINTEND(bugprone-exception-escape)
