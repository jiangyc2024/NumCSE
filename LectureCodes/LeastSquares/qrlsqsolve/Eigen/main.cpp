///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Xiaolin Guo, Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "qrlsqsolve.hpp"
#include <Eigen/Dense>
#include <iostream>


using std::cout;
using std::endl;
using Eigen::Matrix;
using Eigen::VectorXd;

int main()
{
    Matrix<double, 4, 3> A;
    A <<
        1, 2, 3,
        4, 5, 6,
        7, 8, 10,
        1, 1, 1;

    VectorXd b(4);
    b << 3, 3, 4, 0;

    VectorXd x;
    VectorXd x1;
    const double res = qrlsqsolve::qrlsqsolve(A, b, x);
    const double res1 = qrlsqsolve::lsqsolve_eigen(A,b,x1);
    
    cout << "x = " << endl << x << ", res = " << res << endl;
    cout << "x1 = " << endl << x1 << ", res1 = " << res1 << endl;

    return 0;
}
