//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

VectorXd getit(const MatrixXd &A, const VectorXd &x, unsigned int k) {

    EigenSolver<MatrixXd> eig = EigenSolver<MatrixXd>(A);
    const VectorXcd & V = eig.eigenvalues();
    const MatrixXcd & W = eig.eigenvectors();

    VectorXcd cx = x.cast<std::complex<double>>();

    VectorXcd ret = W *
        (
          V.array().pow(k) *
          (W.partialPivLu().solve(cx)).array()
        ).matrix();

    return ret.real();
}

int main() {
    // Some arbitrary data to test getit
    MatrixXd A(4,4);
    A << 1,  2,  3,  4,
         5,  6,  7,  8,
         9,  10, 11, 12,
         13, 14, 15, 16;
    VectorXd x(4);
    x << 4,  5,  6,  7;
    unsigned int  k = 9;

    // Testing the implementation with some matrix
    VectorXd yg = getit(A, x, k);
    std::cout << "getit(A,x, k) = " << std::endl
              << yg << std::endl;

}
