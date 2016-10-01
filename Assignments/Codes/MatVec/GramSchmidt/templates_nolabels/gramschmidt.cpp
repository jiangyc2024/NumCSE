//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

/* \brief Performs Gram-Schidt orthonormalization
 * Given a matrix $\mathbf{A}$ of linearly independent columns,
 * returns the result of a Gram-Schmidt orthonormalization.
 * Unstable GS algorithm: output is prone to cancellation issues.
 * \param[in] $\mathbf{A}$ Matrix of linearly independent columns
 * \return Matrix with ONB of $span(a_1, \cdots, a_n)$ as columns
 */
MatrixXd gram_schmidt(const MatrixXd & A) {
    // We create a matrix Q with the same size and data of A
    MatrixXd Q(A);

    // TODO: Implement GramSchidt othonormalization algorithm

    return Q;
}

int main(void) {
    // Orthonormality test
    unsigned int n = 9;
    MatrixXd A = MatrixXd::Random(n,n);
    // TODO: use gramschmidt to compute orthonormalization of
    // the matrix $\mathbf{A}$.
}
