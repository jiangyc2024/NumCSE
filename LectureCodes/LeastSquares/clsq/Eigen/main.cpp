///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Xiaolin Guo, Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "clsq.hpp"
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;

int main() {
  Eigen::Matrix<double, 4, 4> A;
  // clang-format off
    A <<
      1, 1, 2, 3,
      1, 4, 5, 6,
      1, 7, 8, 10,
      1, 1, 1, 1;
   // clang-format on 
    typename MatrixBase<decltype(A)>::Index dim = 3;
    {
    auto [c,n] = clsq(A, dim);
    cout << "c = " << c << endl;
    cout << "n = " << endl << n << endl;
    }
    {
    auto [c,n] = clsq2(A, 4-dim);
    cout << "c = " << c << endl;
    cout << "n = " << endl << n << endl;
    }
    return 0;
}
