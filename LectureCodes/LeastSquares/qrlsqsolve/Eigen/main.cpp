#include "qrlsqsolve.hpp"
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

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

    VectorXd x,x1;
    double res,res1;
    res = qrlsqsolve(A, b, x);
    res1 = lsqsolve_eigen(A,b,x1);
    
    cout << "x = " << endl << x << ", res = " << res << endl;
    cout << "x1 = " << endl << x1 << ", res1 = " << res << endl;

    return 0;
}
