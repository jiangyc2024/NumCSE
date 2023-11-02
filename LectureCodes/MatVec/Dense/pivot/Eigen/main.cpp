///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>

#include "gausselimsolve.hpp"
#include "lufak.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
  /* SAM_LISTING_BEGIN_0 */
  MatrixXd A(2, 2);
  A << 5.0e-17, 1.0, 1.0, 1.0;
  VectorXd b(2);
  VectorXd x2(2);
  b << 1.0, 2.0;
  const VectorXd x1 = A.fullPivLu().solve(b);
  gausselimsolve::gausselimsolve(A, b, x2); // see Code~\ref{cpp:gausselim}
  const auto [L,U] = lufak::lufak(A);    // see Code~\ref{cpp:lufak}
  const VectorXd z = L.lu().solve(b);
  const VectorXd x3 = U.lu().solve(z);
  std::cout << "x1 = \n"
       << x1 << "\nx2 = \n"
       << x2 << "\nx3 = \n"
       << x3 << std::endl;
  /* SAM_LISTING_END_0 */
  return 0;
}
