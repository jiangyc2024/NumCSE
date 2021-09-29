#ifndef SOLVEPERMB_HPP
#define SOLVEPERMB_HPP

#include <Eigen/Dense>

/**
 * @brief Circular shift (downwards) of $b$
 *
 * @param b The input $n$-dimensional vector shifted downwards
 */
void shift(Eigen::VectorXd& b) {
  const unsigned int n = b.size();

  const double temp = b(n - 1);
  for (int k = n - 2; k >= 0; --k) {
    b(k + 1) = b(k);
  }
  b(0) = temp;
}

/**
 * @brief Compute $X = A^{-1}*[b_1,...,b_n],\; b_i = i$-th cyclic shift of $b$.
 * Function with naive implementation.
 *
 * @param A An $n \times n$ matrix
 * @param b An $n$-dimensional vector
 * @param X The $n \times n$ matrix $X = A^{-1}*[b_1,...,b_n]$
 */
/* SAM_LISTING_BEGIN_0 */
void solvpermb(const Eigen::MatrixXd& A, Eigen::VectorXd& b,
               Eigen::MatrixXd& X) {
  // Size of b, which is the size of A
  const unsigned int n = b.size();
  assert(n == A.rows() && n == A.cols() && "Error: size mismatch!");
  // TODO: (3-11.b) Port the pseudocode from (3.11.1) to C++-code.
  // START

  // END
}
/* SAM_LISTING_END_0 */

/**
 * @brief Compute $X = A^{-1}*[b_1,...,b_n],\; b_i = i$-th cyclic shift of $b$,
 * Function has complexity $O(n^3)$
 *
 * @param A An $n \times n$ matrix
 * @param b An $n$-dimensional vector
 * @param X The $n \times n$ matrix $X = A^{-1}*[b_1,...,b_n]$
 */
/* SAM_LISTING_BEGIN_1 */
void solvpermb_on3(const Eigen::MatrixXd& A, Eigen::VectorXd& b,
                   Eigen::MatrixXd& X) {
  // Size of b, which is the size of A
  const unsigned int n = b.size();
  assert(n == A.cols() && n == A.rows() && "Error: size mismatch!");
  // TODO: (3-11.c) Alter the above function such that it has asymptotic
  // complexity O(n^3).
  // START

  // END
}
/* SAM_LISTING_END_1 */

#endif