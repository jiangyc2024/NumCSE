//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>
#include <cmath>

#include "strassen.hpp"

int main()
{

    // Seed the random number generator
    srand((unsigned int) time(0));

    // Check algorithm for correctness

    // Size of the matrix
    int k = 2;
    int n = std::pow(2, k);

    // Generate random matrices
    MatrixXd A = MatrixXd::Random(n,n);
    MatrixXd B = MatrixXd::Random(n,n);

    // Testing matrix multiplication
    MatrixXd AB = strassenMatMult(A,B);
    // Eigen Mat multiplication
    MatrixXd AxB = A*B;

    std::cout << "Using Strassen's method, A*B=" << std::endl
              << AB << std::endl;
    std::cout << "Using standard method, A*B=" << std::endl
              << AxB << std::endl;
    std::cout << "The norm of the error is: "
              << (AB-AxB).norm() << std::endl;
}
