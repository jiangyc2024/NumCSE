#ifndef KRON_HPP
#define KRON_HPP

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#include "timer.h"

/**
 * @brief Compute the Kronecker product.
 * Computes $\mathbf{C} = \mathbf{A} \otimes \mathbf{B}$.
 * @param A Matrix of size $n \times n$
 * @param B Matrix of size $n \times n$
 * @param C Kronecker product of A and B of dim $n^2 \times n^2$
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd kron(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B) {
  // Allocate enough space for the matrix
  Eigen::MatrixXd C = Eigen::MatrixXd(A.rows() * B.rows(), A.cols() * B.cols());
  // TODO: (1-3.b) Fill in the entries of C.
  // Hint: Use a nested for-loop and C.block().
  // START
  for (unsigned int i = 0; i < A.rows(); ++i) {
    for (unsigned int j = 0; j < A.cols(); ++j) {
      // We use eigen block operations to set the values of
      // each $n \times n$ block.
      C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) =
          A(i, j) * B;  // $\in \mathbb{R}^{(n \times n)}$
    }
  }
  // END

  return C;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Compute the Kronecker product applied to a vector.
 * Computes $\mathbf{y} = (\mathbf{A} \otimes \mathbf{B}) \mathbf{x}$.
 * @param A Matrix of size $n \times n$
 * @param B Matrix of size $n \times n$
 * @param x Vector of dim. $n^2$
 * @param y Vector y = kron(A,B)*x
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd kron_mult(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                          const Eigen::VectorXd &x) {
  assert(A.rows() == A.cols() && A.rows() == B.rows() && B.rows() == B.cols() &&
         "Matrices A and B must be square matrices with same size!");
  assert(x.size() == A.cols() * A.cols() &&
         "Vector x must have length A.cols()^2");
  const unsigned int n = A.rows();

  // Allocate space for output
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n * n);

  // TODO: (1-3.d) Fill in the entries of y.
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
    Eigen::VectorXd z = B * x.segment(j * n, n);
    // Loop over all segments of y
    for (unsigned int i = 0; i < n; ++i) {
      y.segment(i * n, n) += A(i, j) * z;
    }
  }
  // END

  return y;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Compute the Kronecker product $C = A \otimes B$.
 * Use fast reshaping (similar to Matlab reshape)
 * WARNING: using Matrix::Map we assume the matrix is in Column major format,
 *          the code is not valid for Row Major format.
 * @param A Matrix $n \times n$
 * @param B Matrix $n \times n$
 * @param x Vector of dim $n^2$
 * @param y Vector y = kron(A,B)*x
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd kron_reshape(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                             const Eigen::VectorXd &x) {
  assert(A.rows() == A.cols() && A.rows() == B.rows() && B.rows() == B.cols() &&
         "Matrices A and B must be square matrices with same size!");
  const unsigned int n = A.rows();
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n * n);

  // TODO: (1-3.e) Fill in the entires of y.
  // Hint: Use Eigen::MatrixXd::Map() to reshape x into a n by n matrix.
  // Then y is obtained by simple matrix multiplications and
  // another reshape.
  // START
  Eigen::MatrixXd t = B * Eigen::MatrixXd::Map(x.data(), n, n) * A.transpose();
  y = Eigen::MatrixXd::Map(t.data(), n * n, 1);
  // END

  return y;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
void kron_runtime() {
  Eigen::MatrixXd A, B, C;
  Eigen::VectorXd x, y;
  // We repeat each runtime measurement 10 times
  // (this is done in order to remove outliers).
  constexpr unsigned int repeats = 10;

  std::cout << "Runtime for each implementation." << std::endl;
  std::cout << std::setw(5) << "n" << std::setw(15) << "kron" << std::setw(15)
            << "kron_mult" << std::setw(15) << "kron_reshape" << std::endl;
  // Loop from $M = 2,\dots,2^8$
  for (unsigned int M = 2; M <= (1 << 8); M = M << 1) {
    Timer tm_kron, tm_kron_mult, tm_kron_map;
    // Run experiments "repeats" times
    for (unsigned int r = 0; r < repeats; ++r) {
      // Random matrices for testing
      A = Eigen::MatrixXd::Random(M, M);
      B = Eigen::MatrixXd::Random(M, M);
      x = Eigen::VectorXd::Random(M * M);

      // Do not want to use kron for large values of M
      if (M < (1 << 6)) {
        // Kron using direct implementation
        tm_kron.start();
        C = kron(A, B);
        y = C * x;
        tm_kron.stop();
      }

      // TODO: (1-3.f) Measure the runtime of kron_mult() and kron_reshape().
      // START

      // Kron matrix-vector multiplication
      tm_kron_mult.start();
      y = kron_mult(A, B, x);
      tm_kron_mult.stop();

      // Kron using reshape
      tm_kron_map.start();
      y = kron_reshape(A, B, x);
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
