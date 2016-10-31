///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Xiaolin Guo, Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "clsq.hpp"
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main()
{
    Eigen::Matrix<double, 4, 4> A;
    A <<
        1, 1, 2, 3,
        1, 4, 5, 6,
        1, 7, 8, 10,
        1, 1, 1, 1;

    typename MatrixBase<decltype(A)>::Index dim = 3;
    typename MatrixBase<decltype(A)>::Scalar c;
    Eigen::VectorXd n(3);

    clsq(A, dim, c, n);

    cout << "c = " << c << endl;
    cout << "n = " << endl << n << endl;

    return 0;
}