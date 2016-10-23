//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* @brief Build the lower triangular matrix $n \times n$ $A$ from vector $a$.
 * Each column of $A$ contains one element of $a$.
 * \param[in] a An $n$-dimensional vector to build the lower triangular matrix $A$
 * \param[out] A The $n \times n$ lower triangular matrix from $a$
 */
MatrixXd buildA(const VectorXd & a)
{
    // Initialization
    int n = a.size();
    MatrixXd A(n,n);

    for(int j=0; j<n; j++) {
	for(int i=j; i<n; i++) {
	    A(i,j) = a(j);
	}
     }

	return A;
}

/* @brief Solve $Ax = b$.
 * Function with naive implementation.
 * @param[in] a An $n$-dimensional vector to build the lower triangular matrix $A$
 * @param[in] b An $n$-dimensional vector
 * @param[out] x The $n$-dimensional vector $x = A^{-1}*b$
 */
void solveA(const VectorXd & a, const VectorXd & b, VectorXd & x) {
    // Size of b, which is the size of a
    int n = b.size();
    assert( n == a.size()
            && "Error: size mismatch!");
    x.resize(n);

	MatrixXd A = buildA(a);
    x = A.fullPivLu().solve(b);
}

/* @brief Solve $Ax = b$.
 * Function with efficient implementation, using the solution of subproblem 4.
 * @param[in] a An $n$-dimensional vector to build the lower triangular matrix $A$
 * @param[in] b An $n$-dimensional vector
 * @param[out] x The $n$-dimensional vector $x = A^{-1}*b$
 */
void solveA_fast(const VectorXd & a, const VectorXd & b, VectorXd & x) {
    // Size of b, which is the size of a
    int n = b.size();
    assert( n == a.size()
            && "Error: size mismatch!");
    x.resize(n);

    x(0) = b(0)/a(0);
    for(int i=1; i<n; i++) {
		x(i) = (b(i) - b(i-1))/a(i);
	}
}

int main() {
    unsigned int n = 9;
    // Compute with both solvers
    std::cout << "--> Check that the solvers are correct" << std::endl;
    VectorXd a = VectorXd::Random(n);
    VectorXd b = VectorXd::Random(n);
    VectorXd x1, x2, x;

    solveA(a,b,x1);
    solveA_fast(a,b,x2);

    std::cout << "Error = " << (x1 - x2).norm() << std::endl;

    // Compute runtimes of different solvers
    std::cout << "--> Runtime comparison of naive vs fast solver" << std::endl;
	// Number of repetitions
    unsigned int repeats = 3;


    // Header
    std::cout << std::setw(20) << "n"
              << std::setw(20) << "time naive [s]"
              << std::setw(20) << "time fast [s]"
              << std::endl;

    // Loop over matrix size
    for(unsigned int k = 4; k <= 12; ++k) {
        // Timers
        Timer tm_naive, tm_fast;
        unsigned int n = pow(2,k);

        // Repeat test many times
        for(unsigned int r = 0; r < repeats; ++r) {
            a = VectorXd::Random(n);
            b = VectorXd::Random(n);

            // Compute runtime with naive solver
            tm_naive.start();
            solveA(a,b,x);
            tm_naive.stop();
            // Compute runtime with efficient solver
            tm_fast.start();
            solveA_fast(a,b,x);
            tm_fast.stop();
        }


        // Print runtimes
        std::cout << std::setw(20) << n
                  << std::scientific << std::setprecision(3)
                  << std::setw(20) << tm_naive.min()
                  << std::setw(20) << tm_fast.min()
                  << std::endl;
    }

}
