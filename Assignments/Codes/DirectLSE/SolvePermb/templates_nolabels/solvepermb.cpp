//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/LU>

#include "timer.h"

using namespace Eigen;

/* @brief Circular shift (downwards) of b
 * @param[in,out] b The input $n$-dimensional vector shifted downwards
 */
void shift(VectorXd & b) {
    int n = b.size();

    double temp = b(n-1);
    for(int k = n-2; k >= 0; --k) {
        b(k+1) = b(k);
    }
    b(0) = temp;
}

/* @brief Compute $X = inv(A)*[b_1,...,b_n], b_i = i$-th cyclic shift of $b$.
 * Function with naive implementation.
 * @param[in] A An $n \times n$ matrix
 * @param[in] b An $n$-dimensional vector
 * @param[out] X The $n \times n$ matrix $X = inv(A)*[b_1,...,b_n]$
 */
void solvpermb(const MatrixXd & A, VectorXd & b, MatrixXd & X) {
    // Size of b, which is the size of A
    int n = b.size();
    assert( n == A.cols() && n == A.rows()
            && "Error: size mismatch!");
    X.resize(n,n);

    // TODO: solve system $inv(A)*B$
}

/* @brief Compute $X = inv(A)*[b_1,...,b_n], b_i = i$-th cyclic shift of $b$,
 * Function has complexity $O(n^3)$
 * @param[in] A An $n \times n$ matrix
 * @param[in] b An $n$-dimensional vector
 * @param[out] X The $n \times n$ matrix $X = inv(A)*[b_1,...,b_n]$
 */
void solvpermb_on3(const Matrix & A, Vector & b, Matrix & X) {
    // Size of b, which is the size of A
    int n = b.size();
    assert( n == A.cols() && n == A.rows()
            && "Error: size mismatch!");
    X.resize(n,n);

    // TODO: efficiently solve system $inv(A)*B$
}

int main() {
    unsigned int n = 9;
    // Compute with both solvers
    std::cout << "*** Check that the solvers are correct" << std::endl;

    MatrixXd A = MatrixXd::Random(n,n);
    VectorXd b = VectorXd::Random(n);
    MatrixXd X;

    std::cout << "b = " << std::endl
              << b << std::endl;

    solvpermb(A,b,X);
    std::cout << "Direct porting from MATLAB (naive solver): "
              << std::endl << X << std::endl;
    std::cout << "A*X = " << std::endl
              << A*X << std::endl;

    solvpermb_on3(A,b,X);
    std::cout << "Reusing LU: " << std::endl
              << X << std::endl;
    std::cout << "A*X = " << std::endl
              << A*X << std::endl;

    // Compute runtimes of different solvers
    std::cout << "*** Runtime comparison of naive solver vs reusing LU" << std::endl;
    unsigned int repeats = 3;

    // Header
    std::cout << std::setw(10) << "n"
              << std::setw(10) << "time no reuse [s]"
              << std::setw(10) << "time reuse [s]"
              << std::endl;

    for(unsigned int p = 2; p <= 7; ++p) {
        Timer tm_naive, tm_reuseLU;
        unsigned int n = pow(2,p);

        for(unsigned int r = 0; r < repeats; ++r) {
            A = MatrixXd::Random(n,n);
            b = VectorXd::Random(n);

            tm_naive.start();
            solvpermb(A,b,X);
            tm_naive.stop();

            tm_reuseLU.start();
            solvpermb_on3(A,b,X);
            tm_reuseLU.stop();
        }

        std::cout << std::setw(10) << n
                  << std::scientific << std::setprecision(3)
                  << std::setw(10) << tm_naive.min()
                  << std::setw(10) << tm_reuseLU.min()
                  << std::endl;
    }
}
