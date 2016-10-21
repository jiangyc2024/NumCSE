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

/* \brief Performs the computation $y = A^k x$.
 */
VectorXd getit(const MatrixXd &A, const VectorXd &x, unsigned int k) {

    // As said in the problem formulation, we may assume that the
    // following three lines have complexity $O(n^3)$.
    EigenSolver<MatrixXd> eig = EigenSolver<MatrixXd>(A);
    const VectorXcd & V = eig.eigenvalues();
    const MatrixXcd & W = eig.eigenvectors();

    // The following operation requires only a loop over
    // the dimension of $cx$, which is $n$.
    VectorXcd cx = x.cast<std::complex<double>>();

    // The first operator* is a matrix vector multiplication
    // with complexity $O(n^2)$.
    VectorXcd ret = W *
        // The componentwise power has complexity $O(n^2)$.
        // The second operator* is a vector-vector componentwise
        // multiplication, with complexity $O(n)$.
        (
          V.array().pow(k) *
          // In the following line, a linear system is solved, operation
          // with complexity $O(n^3)$
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

    // Checking that getit works
    VectorXd yp = A.pow(k)*x;
    std::cout << "A^k x = " << std::endl
              << yp << std::endl;
    double err = (yg - yp).norm() / yp.norm();
    std::cout << "Relative error = " << err << std::endl;
}
