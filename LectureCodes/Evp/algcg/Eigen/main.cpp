///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting algcg.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>

#include "algcg.hpp"

using namespace Eigen;

int main() {
    // Input values.
    // Assume A is square (so evalA returs a vector size n).
    int n = 2;
    MatrixXd A(n, n);
    std::function<VectorXd(VectorXd)> evalA;
    VectorXd b(n);
    VectorXd x(n);
    double tol = 1e-6;
    unsigned int maxit;

    // Assign arbitrary values.
    A << 1, 3, 5, 7;

    evalA = [A](VectorXd x) { return A * x; };

    b << 1, 2;

    x << 4, 5;

    maxit = 5;

    VectorXd x_approx = cg(evalA, b, x, tol, maxit);

    std::cout << x_approx << "\n";
}
