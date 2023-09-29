///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Xiaolin Guo, Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "lsqsvd.hpp"
#include <Eigen/Dense>
#include <iostream>


using std::cout;
using std::endl;

int main()
{
    Eigen::Matrix<double, 4, 3> A;
    A <<
        1, 2, 3,
        4, 5, 6,
        7, 8, 10,
        1, 1, 1;

    Eigen::VectorXd b(4);
    b << 3, 3, 4, 0;

    cout << lsqsvd::lsqsvd(A, b) << endl;

    return 0;
}