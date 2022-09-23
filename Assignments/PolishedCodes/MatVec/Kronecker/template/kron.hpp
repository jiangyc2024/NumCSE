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

  // TODO: (1-3.e) Fill in the entires of y.
  // Hint: Use Eigen::MatrixXd::Map() to reshape x into a n by n matrix.
  // Then y is obtained by simple matrix multiplications and
  // another reshape.
  // START

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
