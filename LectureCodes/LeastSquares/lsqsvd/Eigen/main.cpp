#include "lsqsvd.hpp"
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

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

    cout << lsqsvd(A, b) << endl;

    return 0;
}