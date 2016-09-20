///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

/* SAM_LISTING_BEGIN_0 */
#include <iostream>
#include <Eigen/Dense>

/**
 *  \brief Given a matrix $A$ of linearly independent columns, returns
 *          Gram-Schmidt orthonormalization
 *
 *  Ustable GS algorithm. Output is prone to cancellation issues.
 *  \param[in] $A$ Matrix of linearly independent columns
 *  \return Matrix with ONB of $span(a_1, \cdots, a_n)$
 */
template <class Matrix>
Matrix gramschmidt( const Matrix & A ) {
    Matrix Q = A;
    // First vector just gets normalized, Line 1 of \eqref{GS}
    Q.col(0).normalize();
    for(unsigned int j = 1; j < A.cols(); ++j) {
        // Replace inner loop over each previous vector in Q with fast
        // matrix-vector multiplication (Lines 4, 5 of \eqref{GS})
        Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));// \Label{gscpp:op}
        // Normalize vector, if possible.
        // Otherwise colums of A must have been linearly dependent
        if( Q.col(j).norm() <= 10e-14 * A.col(j).norm() ) { // \Label{gscpp:1}
            std::cerr << "Gram-Schmidt failed: A has lin. dep columns." << std::endl;
            break;
        } else { Q.col(j).normalize(); } // Line 7 of \eqref{GS}
    }
    return Q;
}
/* SAM_LISTING_END_0 */
