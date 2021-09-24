//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
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
/* SAM_LISTING_BEGIN_0 */
MatrixXd gram_schmidt(const MatrixXd & A) {
    // We create a matrix Q with the same size and data of A
    MatrixXd Q(A);
    
    // The first vector just gets normalized
    Q.col(0).normalize();
    
    // TO DO: (2-2.b) Implement the gram_schmidt procedure by iterating over all other columns of A.
    // START
    for(unsigned int j = 1; j < A.cols(); ++j) {
        // See eigen documentation for usage of col and leftCols
        Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));

        // Normalize vector, if possible
        // (otherwise it means columns of $\mathbf{A}$ are
        // almost linearly dependent)
        double eps = std::numeric_limits<double>::denorm_min();
        if( Q.col(j).norm() <= eps * A.col(j).norm() ) {
            std::cerr << "Gram-Schmidt failed because "
                      << "A has (almost) linearly dependent "
                      << "columns." << std::endl;
            break;
        } else {
            Q.col(j).normalize();
        }
    }
    // END

    return Q;
}
/* SAM_LISTING_END_0 */


/* SAM_LISTING_BEGIN_1 */
double orthogonality_test() {
    // Orthonormality test
    double err=1;
    unsigned int n = 9;
    MatrixXd A, Q;
    A = MatrixXd::Random(n,n);
    
    // TO DO: (2-2.c) Use gram_schmidt() to compute an orthonormalization of A,
    // call it Q, and let err measure "how far Q is from being orthonormal".
    // START
    Q = gram_schmidt( A );
    err = (Q.transpose()*Q - MatrixXd::Identity(n,n)).norm();
    // END
  
    return err;
}
/* SAM_LISTING_END_1 */