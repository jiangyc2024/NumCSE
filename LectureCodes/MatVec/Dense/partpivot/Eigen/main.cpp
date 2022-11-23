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

using namespace std;
using namespace Eigen;

int main() {
  /* SAM_LISTING_BEGIN_0 */
  double epsilon = 5.0e-17;
  MatrixXd A(2, 2), D(2, 2);
  A << epsilon, 1.0, 1.0, 1.0;
  D << 1. / epsilon, 0.0, 0.0, 1.0;
  VectorXd b(2), x2(2);
  b << 1.0, 2.0;
  A = D * A;
  b = D * b;
  VectorXd x1 = A.fullPivLu().solve(b);
  gausselimsolve(A, b, x2); // see Code~\ref{cpp:gausselim}
  auto [L, U] = lufak(A); // see Code~\ref{cpp:lufak}
  VectorXd z = L.lu().solve(b);
  VectorXd x3 = U.lu().solve(z);
  cout << "x1 = \n"
       << x1 << "\nx2 = \n"
       << x2 << "\nx3 = \n"
       << x3 << std::endl;
  /* SAM_LISTING_END_0 */
  return 0;
}
