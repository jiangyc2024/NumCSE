//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>
#include <iomanip>

#include <vector>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* \brief Compute the Kronecker product.
 * Computes $\mathbf{C} = \mathbf{A} \otimes \mathbf{B}$.
 * \param[in] A Matrix of size $n \times n$
 * \param[in] B Matrix of size $n \times n$
 * \param[out] C Kronecker product of A and B of dim $n^2 \times n^2$
 */
void kron(const MatrixXd & A, const MatrixXd & B,
          MatrixXd & C) {
    // Allocate enough space for the matrix
    C = MatrixXd(A.rows()*B.rows(), A.cols()*B.cols());
    for(unsigned int i = 0; i < A.rows(); ++i) {
        for(unsigned int j = 0; j < A.cols(); ++j) {
            // We use eigen block operations to set the values of
            // each $n \times n$ block.
            C.block(i*B.rows(),j*B.cols(), B.rows(), B.cols())
                = A(i,j)*B; // $\in \mathbb{R}^{(n \times n)}$
        }
    }
}

/* \brief Compute the Kronecker product applied to a vector.
 * Computes $\mathbf{y} = (\mathbf{A} \otimes \mathbf{B}) \mathbf{x}$.
 * Exploit efficient matrix-vector product.
 * \param[in] A Matrix of size $n \times n$
 * \param[in] B Matrix of size $n \times n$
 * \param[in] x Vector of dim. $n^2 \times n^2$
 * \param[out] y Vector y = kron(A,B)*x
 */
void kron_mult(const MatrixXd & A, const MatrixXd & B,
               const VectorXd & x, VectorXd & y) {
    assert(A.rows() == A.cols() &&
           A.rows() == B.rows() &&
           B.rows() == B.cols() &&
           "Matrices A and B must be square matrices with same size!");
    assert(x.size() == A.cols()*A.cols() &&
           "Vector x must have length A.cols()^2");
    unsigned int n = A.rows();

    // Allocate space for output
    y = VectorXd::Zero(n*n);


    for(unsigned int j = 0; j < n; ++j) {
        VectorXd z = B * x.segment(j*n, n);
        for(unsigned int i = 0; i < n; ++i) {
            y.segment(i*n, n) += A(i, j)*z;
        }
    }
}

/* \brief Compute the Kronecker product $C = A \otimes B$.
 * Use fast reshaping (similar to Matlab reshape)
 * WARNING: using Matrix::Map we assume the matrix is in Column major format,
 *          the code is not valid for Row Major format.
 * \param[in] A Matrix $n \times n$
 * \param[in] B Matrix $n \times n$
 * \param[in] x Vector of dim $n^2$
 * \param[out] y Vector y = kron(A,B)*x
 */
void kron_reshape(const MatrixXd & A, const MatrixXd & B,
                  const VectorXd & x, VectorXd & y) {
    assert(A.rows() == A.cols() && A.rows() == B.rows() && B.rows() == B.cols() &&
           "Matrices A and B must be square matrices with same size!");
    unsigned int n = A.rows();

    MatrixXd t = B * MatrixXd::Map(x.data(), n, n) * A.transpose();
    y = MatrixXd::Map(t.data(), n*n, 1);
}

int main(void) {

    // Testing correctness of Kron
    MatrixXd A(2,2);
    A << 1, 2, 3, 4;
    MatrixXd B(2,2);
    B << 5, 6, 7, 8;
    MatrixXd C;

    VectorXd x = Eigen::VectorXd::Random(4);
    VectorXd y;

    std::cout << "Testing kron, kron_mult, and kron_reshape with small matrices..."
              << std::endl;

    // Compute using kron
    kron(A, B, C);
    y = C*x;

    std::cout << "kron(A,B) = " << std::endl
              << C << std::endl;
    std::cout << "Using kron: y = " << std::endl
              << y << std::endl;

    // Compute using kron_mult
    kron_mult(A, B, x, y);
    std::cout << "Using kron_mult: y =" << std::endl
              << y << std::endl;

    // Compute using kron_reshape
    kron_reshape(A, B, x, y);
    std::cout << "Using kron_map: y =" << std::endl
              << y << std::endl;

    // Compute runtime of different implementations of kron

    // We repeat each runtime measurment 10 times
    // (this is done in order to remove outliers)
    unsigned int repeats = 10;
    std::vector<double> sizes, times_kron, times_kron_mult, times_kron_map;

    std::cout << "Runtime for each implementation." << std::endl;
    std::cout << std::setw(5) << "n"
              << std::setw(15) << "kron"
              << std::setw(15) << "kron_mult"
              << std::setw(15) << "kron_reshape"
              << std::endl;

    // Loop from $M = 2,\dots,2^5$
    for(unsigned int M = 2; M <= (1 << 8); M = M << 1) {
        Timer tm_kron, tm_kron_mult, tm_kron_map;
        // Run experiments "repeats" times
        for(unsigned int r = 0; r < repeats; ++r) {
            A = MatrixXd::Random(M,M);
            B = MatrixXd::Random(M,M);
            x = VectorXd::Random(M*M);

            if( M < (1 << 6)) {
            tm_kron.start();
                kron(A,B,C);
                y = C*x;
            tm_kron.stop();
            }

            tm_kron_mult.start();
            kron_mult(A,B,x,y);
            tm_kron_mult.stop();

            tm_kron_map.start();
            kron_reshape(A,B,x,y);
            tm_kron_map.stop();
        }

        double kron_time = (M < (1 << 6)) ? tm_kron.min() : std::nan("");
        std::cout << std::setw(5) << M
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << kron_time
                  << std::setw(15) << tm_kron_mult.min()
                  << std::setw(15) << tm_kron_map.min() << std::endl;
        sizes.push_back( M );
        times_kron.push_back( kron_time );
        times_kron_mult.push_back( tm_kron_mult.min() );
        times_kron_map.push_back( tm_kron_map.min() );

    }

}
