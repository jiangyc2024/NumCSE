//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
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
    // TODO: Implement the Kronecker product
    // Hint: Use MatrixXd::block(int, int, int, int)
    // See https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
}

/* \brief Compute the Kronecker product applied to a vector.
 * Computes $\mathbf{y} = (\mathbf{A} \otimes \mathbf{B}) \mathbf{x}$.
 * \param[in] A Matrix of size $n \times n$
 * \param[in] B Matrix of size $n \times n$
 * \param[in] x Vector of dim. $n^2 \times n^2$
 * \param[out] y Vector y = kron(A,B)*x
 */
void kron_mult(const MatrixXd &A, const MatrixXd &B,
               const VectorXd &x, VectorXd &y) {
    assert(A.rows() == A.cols() &&
           A.rows() == B.rows() &&
           B.rows() == B.cols() &&
           "Matrices A and B must be square matrices with same size!");
    assert(x.size() == A.cols()*A.cols() &&
           "Vector x must have length A.cols()^2");
    unsigned int n = A.rows();

    // TODO: implement $y = (A \otimes B) \cdot x$ as efficiently as possible
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

    // TODO: compute y
    // Hint: Use MatrixXd::Map(...)
    // See https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
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
    std::cout << "Using kron_reshape: y =" << std::endl
              << y << std::endl;

    // TODO: Compute runtime of kron, kron_mult and kron_reshape.
    // Output the runtimes in scientific notation with 3 decimal digits
    // Use the Timer class.

}
