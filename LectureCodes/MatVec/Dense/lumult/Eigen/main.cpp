///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>

#include "lumult.hpp"

int main() {
  const int n = 3;
  Eigen::MatrixXd L(n, n);
  Eigen::MatrixXd U(n, n);
  L << 1, 0, 0, -1, 1, 0, -2, 0.25, 1;
  U << 1, 0, 2, 0, 4, 3, 0, 0, 5.25;
  const Eigen::MatrixXd A = lumult::lumult(L, U);

  std::cout << "L=\n" << L << std::endl << "U=\n" << U << std::endl;
  std::cout << "A=\n" << A << std::endl;

  std::cout
      << "|A-L*U| = "
      << (A - (Eigen::MatrixXd(n, n) << 1, 0, 2, -1, 4, 1, -2, 1, 2).finished())
             .norm()
      << std::endl;
  return 0;
}
