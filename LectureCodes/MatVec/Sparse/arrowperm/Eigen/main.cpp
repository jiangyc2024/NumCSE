///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "spy.hpp"
#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::RowVectorXd;

//NOLINTBEGIN(bugprone-exception-escape)
int main() {
  /* SAM_LISTING_BEGIN_0 */
  MatrixXd A(11, 11);
  A.setIdentity();
  A.col(0).setOnes();
  A.row(0) = RowVectorXd::LinSpaced(11, 11, 1);
  // Permutation matrix ($\to$ Def.~\ref{def:permmat}) encoding cyclic
  // permutation
  MatrixXd P(11, 11);
  P.setZero();
  P.topRightCorner(10, 10).setIdentity();
  P(10, 0) = 1;
  spy(A, "A", "InvArrowSpy_cpp.eps");
  spy(P * A * P.transpose(), "permuted A", "ArrowSpy_cpp.eps");
  /* SAM_LISTING_END_0 */
  return 0;
}
//NOLINTEND(bugprone-exception-escape)