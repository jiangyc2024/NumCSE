#include <iostream>

/**
 *  \brief Given a matrix $A$ of linearly independent columns, returns 
 *          Gram-Schmidt orthonormalization
 *   
 *  Ustable GS algorithm. Output is prone to cancellation issues.
 *  \param[in] $A$ Matrix of linearly independent columns
 *  \return Matrix with ONB of $span(a_1, \cdots, a_n)$
 */
template <class Matrix>
Matrix gramschmidt( const Matrix & A ) { /* SAM_SOLUTION_BEGIN */
    Matrix Q(A);

    // First vector just gets normalized
    Q.col(0).normalize();
    
    for(unsigned int j = 1; j < A.cols(); ++j) {
        // Replace inner loop over each previous vector in Q with fast 
        //  matrix-vector multiplication
        Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));
        
        // Normalize vector if possible (otherwise means colums of $A$
        //  almost linear dependant)
        if( Q.col(j).norm() <= 10e-14 * A.col(j).norm() ) {
            std::cerr << "Gram-Schmidt failed because A has linear dependant"
                      << " columns. Bye." << std::endl;
            break;
        } else {
            Q.col(j).normalize();
        }
    }
    
    return Q; /* SAM_SOLUTION_END */
}
