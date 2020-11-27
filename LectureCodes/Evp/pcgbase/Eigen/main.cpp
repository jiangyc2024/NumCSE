///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting pcgbase.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "pcgbase.hpp"

int main() {
    // NOTE: C++ implementation returns quite different values than MATLAB. 
    //       Likely rounding/cancellation errors.
    // Declare variables.
    // Assume A is square.
    int n = 2;
    MatrixXd A(n, n);
    std::function<VectorXd(VectorXd)> evalA;
    VectorXd b(n);
    double tol = 1e-6;
    unsigned int maxit;
    std::function<VectorXd(VectorXd)> invB;
    VectorXd x(n);

    // Initialise variables arbitrarily.
    A << 1, 2,
         3, 4;

    evalA = [A](VectorXd x) {
        return A * x;
    };

    b << 1, 2;
    maxit = 5;

    invB = [](VectorXd x) {
        return x;
    };

    x << 1, 2;

    // Output
    auto out = pcgbase(evalA, b, tol, maxit, invB, x);
    auto x_out = std::get<0>(out);
    auto rn_out = std::get<1>(out);
    auto xk_out = std::get<2>(out);

    std::cout << "x_out:\n"
              << "(" << x_out.transpose() << ")";
    std::cout << "\n\nrn_out:\n";
    for (int i = 0; i < rn_out.size(); i++) {
        std::cout << rn_out[i] << "\n";
    }
    std::cout << "\n\nxk_out:\n";
    for (int i = 0; i < xk_out.size(); i++) {
        std::cout << "(" << xk_out[i].transpose() << ")\n";
    }
}