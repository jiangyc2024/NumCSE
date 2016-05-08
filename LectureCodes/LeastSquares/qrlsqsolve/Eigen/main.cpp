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

    VectorXd x;
    double res;
    qrlsqsolve(A, b, x, res);

    cout << "x = " << endl << x << endl;
    cout << "res = " << res << endl;

    return 0;
}