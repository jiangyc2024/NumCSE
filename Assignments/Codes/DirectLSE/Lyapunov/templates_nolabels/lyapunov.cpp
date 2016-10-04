//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

/* @brief Compute the matrix $C$ from $A$
 * @param[in] A An $n \times n$ matrix
 * @param[out] C The $(n^2) \times (n^2)$ matrix
 * from $A\otimes I+I\otimes A$
 */
SparseMatrix<double> buildC(const MatrixXd &A)
{
    // Initialization
    int n = A.rows();
    SparseMatrix<double> C(n*n,n*n);
    std::vector<Triplet<double> > triplets;
    MatrixXd I = MatrixXd::Identity(n,n);

    // TODO: compute $C$ from $A$

	return C;
}

/* @brief Solve the Lyapunov system
 * @param[in] A An $n \times n$ matrix
 * @param[out] X The $n \times n$ solution matrix
 */
void solveLyapunov(const MatrixXd &A, MatrixXd &X)
{
    // Initialization
    int n = A.rows();
    SparseMatrix <double> C = buildC(A);
    MatrixXd I = MatrixXd::Identity(n,n);
    VectorXd b(n*n);
    b = Map<MatrixXd>(I.data(),n*n,1);
    VectorXd vecX(n*n);

    // TODO: solve Lyapunov
}

int main() {
    // Initialization
    unsigned int n = 5;
    MatrixXd A(n,n);
    A << 10, 2, 3, 4, 5,
         6, 20, 8, 9, 1,
         1, 2, 30, 4, 5,
         6, 7, 8, 20, 0,
         1, 2, 3, 4, 10;

    // Test 'buildC'
    SparseMatrix <double> C = buildC(A);
    std::cout << "C = " << C << std::endl;

    // Test 'solveLyapunov'
    MatrixXd X(n,n);
    solveLyapunov(A,X);
    std::cout << "X = " << X << std::endl;

    // Verify the solution if you obtain zero
    MatrixXd I = MatrixXd::Identity(n,n);
    std::cout << "Correct if close to 0: "
              << (A*X + X*A.transpose() - I).norm()
              << std::endl;
}
