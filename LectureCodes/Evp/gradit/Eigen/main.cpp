///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting gradit.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////


#include "gradit.hpp"
#include <Eigen/Dense>
#include <iostream>

int main() {
    // Input values. 
    // b and x must be column vectors of same size.
    // Asume A is square.

    int n = 2;
    MatrixXd A(n, n);
    VectorXd b(n);
    VectorXd x(n);
    double rtol = 1e-6;
    double atol = 1e-6;
    unsigned int maxit = 5;

    // Assign random values.
    A << 1, 2,
         3, 4;
    b << 10, 11;
    x << 20, 21;

    // Output
    VectorXd x_it = gradit(A, b, x, rtol, atol, maxit);
    std::cout << x_it << "\n";
}