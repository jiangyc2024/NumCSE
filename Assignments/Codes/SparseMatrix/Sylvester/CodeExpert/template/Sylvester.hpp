#include <iomanip>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

/* @brief Solve $X+AXA=A$ for $X$ where $A$ is a diagonal s.p.d. matrix.
 * \param[in] diagA, a size $n$ Vector of the diagonal entries of A (all > 0).
 * \param[out] X, the $n\times n$ sparse matrix that solves $X+AXA=A$.
 */
/* SAM_LISTING_BEGIN_0 */
template <typename Vector>
SparseMatrix<double> solveDiagSylvesterEq(const Vector &diagA) {
  int n = diagA.size();
  SparseMatrix<double> X(n, n);
  // TO DO: Fill in the entries of X.
  // Don't forget to use makeCompressed().
  // START
  
  // END
  return X;
}
/* SAM_LISTING_END_0 */

/* @brief Compute the Kronecker product $A\otimes A$.
 * \param[in] A, an $n\times n$ sparse matrix.
 * \param[out] B, the $n^2\times n^2$ sparse matrix $B=A\otimes A$
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> sparseKron(const SparseMatrix<double> &M) {
  int n = M.rows();
  assert(n == M.cols() && "Matrix M must be square.");

  // JM is an array of size n+1, with
  // JM[j] = number of non-zero entries in all columns 0,..., (j-1) of M.
  const int *JM = M.outerIndexPtr();

  // Define nnzColwiseM(j) = nnz( M.col(j) )
  VectorXd nnzColwiseM(n);
  for (int j = 0; j < n; j++) {
    nnzColwiseM(j) = JM[j + 1] - JM[j];
  }

  // B will be the Kronecker product of M with itself.
  SparseMatrix<double> B(n * n, n * n);

  // TO DO: Fill in the entries of B.
  // Hint: You can use nnzColwiseM when reserving space for B.
  // Use M.valuePtr() and M.innerIndexPtr() to define arrays similar
  // to JM. Use those arrays to access the non-zero entries of M.
  // START
  
  // END
  return B;
}
/* SAM_LISTING_END_1 */

/* @brief Solve $XA^{-1}+AX=I$ for $X$ where $A$ is s.p.d.
 * \param[in] A, sparse s.p.d. matrix.
 * \param[out] X, the $n\times n$ sparse matrix that solves $X+AXA=A$.
 */
/* SAM_LISTING_BEGIN_2 */
MatrixXd solveSpecialSylvesterEq(const SparseMatrix<double> &A) {
  const int n = A.rows();
  assert(n == A.cols() && "Matrix A must be square.");
  MatrixXd X(n, n);

  // TO DO: Solve the equation X*A^{-1} + A*X = I .
  // START
  
  // END
  return X;
}
/* SAM_LISTING_END_2 */
