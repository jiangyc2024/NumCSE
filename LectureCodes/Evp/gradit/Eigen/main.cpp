///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting gradit.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson, Vivienne Langen
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

    // !!  OPEN QUESTION !! 
    // regarding: "..should just be operator(), no more evalA.."
    auto evalA = [A](VectorXd x) { return A * x; };


    // Output
    VectorXd x_it = gradit(A, b, x, rtol, atol, maxit);
    std::cout << x_it << "\n";
}
