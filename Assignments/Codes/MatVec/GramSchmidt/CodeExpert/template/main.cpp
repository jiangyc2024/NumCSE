//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>

#include "gramschmidt.hpp"

using namespace Eigen;


int main(void) {
    double tolerance = 1e-9;
    Matrix4d A;
    A << 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1;
    std::cout << "Is correct for test matrix: "
            << ((gram_schmidt(A) - Matrix4d::Identity()).norm() < tolerance)
            << std::endl;
  
    // Orthonormality test
    double err;
    std::srand(5);
    err = orthogonality_test();
    
    // Error has to be small, but not zero (why?)
    std::cout << "Error is: "
              << err
              << std::endl;
    
    std::cout << "Is orthogonal: " << ((err < tolerance) && (err >= 0)) << std::endl;


}
