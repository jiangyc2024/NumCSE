///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Xiaolin Guo, Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "gn.hpp"
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main()
{
    //////////////////////////////////////////
    // solve linear least square problem Ax=b
    //////////////////////////////////////////
    cout << "solving Ax=b" << endl;

    Eigen::Matrix<double, 4, 3> A;
    A <<
        1, 2, 3,
        4, 5, 6,
        7, 8, 10,
        1, 1, 1;
    Eigen::VectorXd b(4);
    b << 3, 3, 4, 0;
    Eigen::VectorXd x0(3);
    x0 << 1, 1, 1;

    // argmin ||F(x)|| = argmin ||Ax-b||
    function<VectorXd(const VectorXd &)> F1 =
        [&A, &b](VectorXd x)
        { return A * x - b; };
    // the Jacobian matrix of F at x is A
    function<MatrixXd(const VectorXd &)> J1 =
        [&A](VectorXd x)
        { return A; };

    cout << gn(x0, F1, J1, 0.001) << endl;

    //////////////////////////////////////////
    // solve nonlinear least square problem
    // see https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm#Example
    //////////////////////////////////////////
    cout << "solving nonlinear least square" << endl;

    VectorXd X(7), Y(7);
    X << 0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740;
    Y << 0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317;

    function<VectorXd(const VectorXd &)> F2 =
        [&X, &Y](VectorXd beta)
        { return beta(0) * X.array() / (beta(1) + X.array()) - Y.array(); };

    function<MatrixXd(const VectorXd &)> J2 =
        [&X](VectorXd beta)
        {
            MatrixXd J(X.rows(), beta.rows());
            J.col(0) = X.array() / (beta(1) + X.array());
            J.col(1) = -beta(0) * X.array() / (beta(1) + X.array()).square();
            return J;
        };

    VectorXd beta0(2);
    beta0 << 0.9, 0.2;
    cout << gn(beta0, F2, J2, 1e-6) << endl;

    return 0;
}
