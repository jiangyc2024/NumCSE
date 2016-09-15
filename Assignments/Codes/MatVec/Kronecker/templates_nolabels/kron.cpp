//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "timer.hpp"

using Eigen;

/* \brief Compute the Kronecker product.
 * Computes $\mathbf{C} = \mathbf{A} \otimes \mathbf{B}$.
 * \param[in] A Matrix of size $n \times n$
 * \param[in] B Matrix of size $n \times n$
 * \param[out] C Kronecker product of A and B of dim $n^2 \times n^2$
 */
void kron(const MatrixXd & A, const MatrixXd & B,
          MatrixXd & C) {
    // TODO: Implement the Kronecher product
    // Hint: Use MatrixXd::block(int, int, int, int)
    // See https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
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
    assert(A.rows() == A.cols() && A.rows() == B.rows() && B.rows() == B.cols() &&
           "Matrices A and B must be square matrices with same size!").
    unsigned int n = A.rows();

    // TODO: implement $y = (A \otimes B) \cdot x$
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
           "Matrices A and B must be square matrices with same size!").
    unsigned int n = A.rows();

    // TODO: compute y using MatrixXd::Map();
    // Hint: Use MatrixXd::Map(...)
    // See https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
}

int main(void) {

    // Check if kron works, cf.
    Eigen::MatrixXd A(2,2);
    A << 1, 2, 3, 4;
    Eigen::MatrixXd B(2,2);
    B << 5, 6, 7, 8;
    Eigen::MatrixXd C;

    Eigen::VectorXd x = Eigen::VectorXd::Random(4);
    Eigen::VectorXd y;
    kron(A,B,C);
    y = C*x;
    std::cout << "kron(A,B)=" << std::endl << C << std::endl;
    std::cout << "Using kron: y=       " << std::endl << y << std::endl;

    kron_fast(A,B,x,y);
    std::cout << "Using kron_fast: y=  " << std::endl << y << std::endl;
    kron_super_fast(A,B,x,y);
    std::cout << "Using kron_super_fast: y=  " << std::endl << y << std::endl;

    // Compute runtime of different implementations of kron
    unsigned int repeats = 10;
    std::vector<double> times_kron, times_kron_fast, times_kron_super_fast;

    for(unsigned int p = 2; p <= 9; p++) {
        Timer tm_kron, tm_kron_fast, tm_kron_super_fast;
        for(unsigned int r = 0; r < repeats; ++r) {
            unsigned int M = pow(2,p);
            A = Eigen::MatrixXd::Random(M,M);
            B = Eigen::MatrixXd::Random(M,M);
            x = Eigen::VectorXd::Random(M*M);

            // May be too slow for large p, comment if so
            tm_kron.start();
            //     kron(A,B,C);
            //     y = C*x;
            tm_kron.stop();

            tm_kron_fast.start();
            kron_fast(A,B,x,y);
            tm_kron_fast.stop();

            tm_kron_super_fast.start();
            kron_super_fast(A,B,x,y);
            tm_kron_super_fast.stop();
        }

        std::cout << "Lazy Kron took:       " << tm_kron.min() << " s" << std::endl;
        std::cout << "Kron fast took:       " << tm_kron_fast.min() << " s" << std::endl;
        std::cout << "Kron super fast took: " << tm_kron_super_fast.min() << " s" << std::endl;
        times_kron.push_back( tm_kron.min() );
        times_kron_fast.push_back( tm_kron_fast.min() );
        times_kron_super_fast.push_back( tm_kron_super_fast.min() );
    }

    for(auto it = times_kron.begin(); it != times_kron.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    for(auto it = times_kron_fast.begin(); it != times_kron_fast.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    for(auto it = times_kron_super_fast.begin(); it != times_kron_super_fast.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
}
