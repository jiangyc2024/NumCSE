#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "timer.h"

//! \brief Given a matrix A of linearly independent columns, returns Gram-Schmidt orthonormalization
//! \param[in] A Matrix of linearly independent columns
//! \return Matrix with ONB of $span(a_1, \cdots, a_n)$
template <class Matrix>
Matrix gramschmidt( const Matrix & A ) {
    
    Matrix Q = A;
    Q.col(0).normalize();
    
    for(unsigned int j = 1; j < A.cols(); ++j) {
        for(unsigned int l = 0; l < j; ++l) {
            Q.col(j) -= A.col(j).dot(Q.col(l)) * Q.col(l);
        }
        if( Q.col(j).norm() <= 10e-14 * A.col(j).norm() ) {
            std::cerr << "Gram-Schmidt failed because A has lin. dep columns. Bye." << std::endl;
            break;
        } else {
            Q.col(j).normalize();
        }
    }
    
    return Q;
}

//! \brief Given a matrix A of linearly independent columns, returns Gram-Schmidt orthonormalization
//! Version using no internal loop, instead uses matrix-matrix mult.
//! \param[in] A Matrix of linearly independent columns
//! \return Matrix with ONB of $span(a_1, \cdots, a_n)$
template <class Matrix>
Matrix gramschmidt_lessloops( const Matrix & A ) {
    
    Matrix Q = A;
    Q.col(0).normalize();
    
    for(unsigned int j = 1; j < A.cols(); ++j) {
        Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));
        
        if( Q.col(j).norm() <= 10e-14 * A.col(j).norm() ) {
            std::cerr << "Gram-Schmidt failed because A has lin. dep columns. Bye." << std::endl;
            break;
        } else {
            Q.col(j).normalize();
        }
    }
    
    return Q;
}

//! \brief Generates Hilbert matrix, given a matrix H with defined shape
//! \param[in,out] Input matrix for the shape, output written into H
template <class Matrix>
void hilbert( Matrix & H ) {
    for(unsigned i = 0; i < H.rows(); ++i) {
        for(unsigned j = 0; j < H.cols(); ++j) {
            H(i,j) = 1. / (i+j+1);
        }
    }
}

int main(void) {
    // Ortho test
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3,3);
    Eigen::MatrixXd Q = gramschmidt( A );
//     Eigen::MatrixXd Q = gramschmidt_lessloops( A );
    
    // Output should be idenity matrix
    std::cout << Q*Q.transpose() << std::endl;
    
    timer<> tm, tm_fast;
    std::vector<int> times, times_fast;
    unsigned N = 9, repeats = 3, p = 1;
    for(unsigned n = 0; n < N; ++n) {
        for(unsigned int r = 0; r < repeats; ++r) {
            A = Eigen::MatrixXd::Random(p,p);
            tm.start();
            Q = gramschmidt( A );
            tm.stop();
            tm_fast.start();
            Q = gramschmidt_lessloops( A );
            tm_fast.stop();
        }
        times.push_back( tm.avg().count() );
        times_fast.push_back( tm_fast.avg().count() );
        p *= 2;
    }
    for(auto it = times.begin(); it != times.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    for(auto it = times_fast.begin(); it != times_fast.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    
    A = Eigen::MatrixXd::Zero(10,10);
    
    hilbert(A);
    std::cout << A << std::endl;
    Q = gramschmidt( A );
    std::cout << Q*Q.transpose() << std::endl;
    Q = gramschmidt_lessloops( A );
    std::cout << Q*Q.transpose() << std::endl;
}
