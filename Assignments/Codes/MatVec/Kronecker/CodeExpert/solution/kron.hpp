#ifndef KRON_HPP
#define KRON_HPP


#include <iomanip>
#include <iostream>

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
/* SAM_LISTING_BEGIN_1 */
void kron(const MatrixXd &A, const MatrixXd &B, MatrixXd &C) {
  // Allocate enough space for the matrix
  C = MatrixXd(A.rows() * B.rows(), A.cols() * B.cols());
  // TO DO: (2-3.b) Fill in the entries of C.
  // Hint: Use a nested for-loop and C.block().
  // START
  for (unsigned int i = 0; i < A.rows(); ++i) {
    for (unsigned int j = 0; j < A.cols(); ++j) {
      // We use eigen block operations to set the values of
      // each $n \times n$ block.
      C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) =
          A(i, j) * B; // $\in \mathbb{R}^{(n \times n)}$
    }
  }
  // END
}
/* SAM_LISTING_END_1 */

/* \brief Compute the Kronecker product applied to a vector.
 * Computes $\mathbf{y} = (\mathbf{A} \otimes \mathbf{B}) \mathbf{x}$.
 * \param[in] A Matrix of size $n \times n$
 * \param[in] B Matrix of size $n \times n$
 * \param[in] x Vector of dim. $n^2 \times n^2$
 * \param[out] y Vector y = kron(A,B)*x
 */
/* SAM_LISTING_BEGIN_2 */
void kron_mult(const MatrixXd &A, const MatrixXd &B, const VectorXd &x,
               VectorXd &y) {
  assert(A.rows() == A.cols() && A.rows() == B.rows() && B.rows() == B.cols() &&
         "Matrices A and B must be square matrices with same size!");
  assert(x.size() == A.cols() * A.cols() &&
         "Vector x must have length A.cols()^2");
  unsigned int n = A.rows();

  // Allocate space for output
  y = VectorXd::Zero(n * n);
  
  // TO DO: (2-3.d) Fill in the entires of y.
  // Hint: Use a nested for-loop, x.segment(), and y.segment().
  // In the outer loop, you can perform a computation based on
  // B and x, and save the result in a variable that is reused in
  // each iteration of the inner loop.
  // START
  
  // Note: this is like doing a matrix-vector multiplication
  // where the entries of the matrix are smaller matrices
  // and entries of the vector are smaller vectors

  // Loop over all segments of x ($\tilde{x}$)
  for (unsigned int j = 0; j < n; ++j) {
    // Reuse computation of z
    VectorXd z = B * x.segment(j * n, n);
    // Loop over all segments of y
    for (unsigned int i = 0; i < n; ++i) {
      y.segment(i * n, n) += A(i, j) * z;
    }
  }
  // END
}
/* SAM_LISTING_END_2 */

/* \brief Compute the Kronecker product $C = A \otimes B$.
 * Use fast reshaping (similar to Matlab reshape)
 * WARNING: using Matrix::Map we assume the matrix is in Column major format,
 *          the code is not valid for Row Major format.
 * \param[in] A Matrix $n \times n$
 * \param[in] B Matrix $n \times n$
 * \param[in] x Vector of dim $n^2$
 * \param[out] y Vector y = kron(A,B)*x
 */
/* SAM_LISTING_BEGIN_3 */
void kron_reshape(const MatrixXd &A, const MatrixXd &B, const VectorXd &x,
                  VectorXd &y) {
  assert(A.rows() == A.cols() && A.rows() == B.rows() && B.rows() == B.cols() &&
         "Matrices A and B must be square matrices with same size!");
  unsigned int n = A.rows();
  
  // TO DO: (2-3.e) Fill in the entires of y.
  // Hint: Use MatrixXd::Map() to reshape x into a n by n matrix.
  // Then y is obtained by simple matrix multiplications and
  // another reshape.
  // START
  MatrixXd t = B * MatrixXd::Map(x.data(), n, n) * A.transpose();
  y = MatrixXd::Map(t.data(), n * n, 1);
  // END
}
/* SAM_LISTING_END_3 */


/* SAM_LISTING_BEGIN_4 */
void kron_runtime() {
  
  MatrixXd A, B, C;
  VectorXd x, y;
  // We repeat each runtime measurement 10 times
  // (this is done in order to remove outliers).
  unsigned int repeats = 10;

  std::cout << "Runtime for each implementation." << std::endl;
  std::cout << std::setw(5) << "n" << std::setw(15) << "kron" << std::setw(15)
            << "kron_mult" << std::setw(15) << "kron_reshape" << std::endl;
  // Loop from $M = 2,\dots,2^8$
  for (unsigned int M = 2; M <= (1 << 8); M = M << 1) {
    Timer tm_kron, tm_kron_mult, tm_kron_map;
    // Run experiments "repeats" times
    for (unsigned int r = 0; r < repeats; ++r) {
      // Random matrices for testing
      A = MatrixXd::Random(M, M);
      B = MatrixXd::Random(M, M);
      x = VectorXd::Random(M * M);
      
      
      // Do not want to use kron for large values of M
      if (M < (1 << 6)) {
        // Kron using direct implementation
        tm_kron.start();
        kron(A, B, C);
        y = C * x;
        tm_kron.stop();
      }
      
      // TO DO: (2-3.f) Measure the runtime of kron_mult()
      // and kron_reshape().
      // START

      // Kron matrix-vector multiplication
      tm_kron_mult.start();
      kron_mult(A, B, x, y);
      tm_kron_mult.stop();

      // Kron using reshape
      tm_kron_map.start();
      kron_reshape(A, B, x, y);
      tm_kron_map.stop();
      
      // END
    }

    double kron_time = (M < (1 << 6)) ? tm_kron.min() : std::nan("");
    std::cout << std::setw(5) << M << std::scientific << std::setprecision(3)
              << std::setw(15) << kron_time << std::setw(15)
              << tm_kron_mult.min() << std::setw(15) << tm_kron_map.min()
              << std::endl;
  }
}
/* SAM_LISTING_END_4 */

#endif