#ifndef STRUCTUREDLSE_HPP
#define STRUCTUREDLSE_HPP

#include <Eigen/Dense>
#include <cassert>

/**
 * @brief Build the lower triangular matrix $n \times n$ $A$ from vector $a$.
 * Each column of $A$ contains one element of $a$.
 *
 * @param a An $n$-dimensional vector to build the lower triangular matrix $A$
 * @return Eigen::MatrixXd The $n \times n$ lower triangular matrix from $a$
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd buildA(const Eigen::VectorXd& a) {
  // Initialization
  const unsigned int n = a.size();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
  // TODO: (3-12.b) Build the matrix A given the vector a.
  // START
  for (unsigned int j = 0; j < n; ++j) {
    for (unsigned int i = j; i < n; ++i) {
      A(i, j) = a(j);
    }
  }
  // END
  return A;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Solve $Ax = b$.
 * Function with naive implementation.
 *
 * @param a An $n$-dimensional vector to build the lower triangular matrix $A$
 * @param b An $n$-dimensional vector
 * @param x The $n$-dimensional vector $x = A^{-1}*b$
 */
/* SAM_LISTING_BEGIN_1 */
void solveA(const Eigen::VectorXd& a, const Eigen::VectorXd& b,
            Eigen::VectorXd& x) {
  // size of b, which is the size of a
  const unsigned int n = b.size();
  assert(n == a.size() && "Error: size mismatch!");
  // TODO: (3-12.c) Solve Ax = b using "structure oblivious" Gaussian
  // elimination.
  // START
  x.resize(n);
  Eigen::MatrixXd A = buildA(a);
  x = A.fullPivLu().solve(b);
  // END
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solve $Ax = b$.
 * Function with efficient implementation, using the solution of subproblem
 * 3-12.d.
 *
 * @param a An $n$-dimensional vector to build the lower triangular matrix $A$
 * @param b An $n$-dimensional vector
 * @param x The $n$-dimensional vector $x = A^{-1}*b$
 */
/* SAM_LISTING_BEGIN_2 */
void solveA_fast(const Eigen::VectorXd& a, const Eigen::VectorXd& b,
                 Eigen::VectorXd& x) {
  // size of b, which is the size of a
  const unsigned int n = b.size();
  assert(n == a.size() && "Error: size mismatch!");
  // TODO: (3-12.e) Efficiently (O(n)) solve Ax = b.
  // START
  x.resize(n);

  x(0) = b(0) / a(0);
  for (unsigned int i = 1; i < n; ++i) {
    x(i) = (b(i) - b(i - 1)) / a(i);
  }
  // END
}
/* SAM_LISTING_END_2 */

#endif