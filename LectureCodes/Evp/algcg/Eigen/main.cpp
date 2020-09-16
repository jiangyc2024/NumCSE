///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting algcg.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////


#include "algcg.hpp"
#include <Eigen/Dense>
#include <iostream>

int main() {
    // Input values.
    // Assume A is square (so evalA returs a vector size n).
    int n = 2;
    // TODO: evalA
    VectorXd b(n);
    VectorXd x(n);
    double tol = 1e-6;
    unsigned int maxit;

    // Assign arbitrary values.
    // TODO: evalA
    b << 1,
         2;

    x << 4,
         5;
    
    maxit = 5;

    VectorXd x_approx = cg(evalA, b, x, tol, maxit);

    std::cout << x_approx << "\n";
}