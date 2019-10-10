#ifndef LYAPUNOV_HPP
#define LYAPUNOV_HPP

#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

/* @brief Compute the matrix $C$ from $A$
 * @param[in] A An $n \times n$ matrix
 * @param[out] C The $(n^2) \times (n^2)$ matrix
 * from $A\otimes I+I\otimes A$
 */
/* SAM_LISTING_BEGIN_0 */
SparseMatrix<double> buildC(const MatrixXd &A) {
  // Initialization
  int n = A.rows();
  SparseMatrix<double> C(n * n, n * n);
  std::vector<Triplet<double>> triplets;
  MatrixXd I = MatrixXd::Identity(n, n);

  // Iterate over n*n blocks.
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // TO DO: (3-8.e) Fill in the (i,j)-th block of C.
      // Hint: The vector triplets should contain the position
      // and value the non-zero entries of C.
      // Use triplets.push_back() to add a new
      // Triplet<double> to the list.
      // START
      
      // END
    }
  }
  C.setFromTriplets(triplets.begin(), triplets.end());
  C.makeCompressed();
  return C;
}
/* SAM_LISTING_END_0 */

/* @brief Solve the Lyapunov system
 * @param[in] A An $n \times n$ matrix
 * @param[out] X The $n \times n$ solution matrix
 */
/* SAM_LISTING_BEGIN_1 */
void solveLyapunov(const MatrixXd &A, MatrixXd &X) {
  // To solve A*X + X*A^T = I, we must first write
  // the LSE in "canonical form", C*vecX = b.
  int n = A.rows();
  SparseMatrix<double> C;
  MatrixXd I = MatrixXd::Identity(n, n);
  VectorXd b(n * n);
  VectorXd vecX(n * n);

  // TO DO: (3-8.f) Fill in the entries of C, b, and vecX.
  // Hint: Use SparseLU to solve the LSE.
  // START
  
  // END

  X = Map<MatrixXd>(vecX.data(), n, n);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
bool testLyapunov() {
  // Initialization
  unsigned int n = 5;
  MatrixXd A(n, n);
  A << 10, 2, 3, 4, 5, 6, 20, 8, 9, 1, 1, 2, 30, 4, 5, 6, 7, 8, 20, 0, 1, 2, 3,
      4, 10;

  double tol = 1e-8;
  bool works;

  // TO DO: (3-8.g)
  // i) Print the full system matrix C (with zeroes),
  // ii) print the solution X, and
  // iii) return true if and only if the error defined
  // by || A*X + X*A^T - I || is less than tol.
  // START
  
  // END

  return works;
}
/* SAM_LISTING_END_2 */

#endif
