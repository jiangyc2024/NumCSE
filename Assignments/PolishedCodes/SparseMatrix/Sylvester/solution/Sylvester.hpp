#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>

/* @brief Solve $X+AXA=A$ for $X$ where $A$ is a diagonal s.p.d. matrix.
 * @param diagA, a size $n$ Vector of the diagonal entries of A (all > 0).
 * @param X, the $n\times n$ sparse matrix that solves $X+AXA=A$.
 */
/* SAM_LISTING_BEGIN_0 */
template <typename Vector>
Eigen::SparseMatrix<double> solveDiagSylvesterEq(const Vector &diagA) {
  int n = diagA.size();
  Eigen::SparseMatrix<double> X(n, n);

  // TO DO: Fill in the entries of X.
  // Don't forget to use makeCompressed().

  // START
  // As X is a diagonal matrix, we reserve space
  // for one non-zero entry per column.
  X.reserve(Eigen::VectorXi::Constant(n, 1));
  for (int j = 0; j < n; j++) {
    // Set diagonal entry $\cob{\VX_{j,j} = \frac{d_j}{d_j^2 + 1}}$
    X.insert(j, j) = diagA[j] / (diagA[j] * diagA[j] + 1.0);
  }
  X.makeCompressed();
  // END

  return X;
}
/* SAM_LISTING_END_0 */

/* @brief Compute the Kronecker product $A\otimes A$.
 * @param A, an $n\times n$ sparse matrix.
 * @param B, the $n^2\times n^2$ sparse matrix $B=A\otimes A$
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> sparseKron(const Eigen::SparseMatrix<double> &M) {
  int n = M.rows();
  assert(n == M.cols() && "Matrix M must be square.");

  // JM is an array of size n+1, with
  // JM[j] = number of non-zero entries in all columns 0,..., (j-1) of M.
  const int *JM = M.outerIndexPtr();

  // Define nnzColwiseM(j) = nnz( M.col(j) )
  Eigen::VectorXi nnzColwiseM(n);
  for (int j = 0; j < n; j++) {
    nnzColwiseM(j) = JM[j + 1] - JM[j];
  }

  // B will be the Kronecker product of M with itself.
  Eigen::SparseMatrix<double> B(n * n, n * n);

  // TO DO: Fill in the entries of B.
  // Hint: You can use nnzColwiseM when reserving space for B.
  // Use M.valuePtr() and M.innerIndexPtr() to define arrays similar
  // to JM. Use those arrays to access the non-zero entries of M.

  // START

  Eigen::VectorXi nnzColwiseB;
  // We have that nnz(B.col(j*n+l)) = nnz(M.col(j)) * nnz(M.col(l)).
  Eigen::MatrixXi tmp = nnzColwiseM * nnzColwiseM.transpose();  // O(n^2)
  nnzColwiseB = Eigen::MatrixXi::Map(tmp.data(), n * n, 1);

  // Reserve sufficient space.
  B.reserve(nnzColwiseB);

  // We have that B.block(i*n,j*n,n,n) = M(i,j) * M,
  // i.e. B( i*n + k, j*n + l ) = M(i,j)*M(k,l).
  // Hence, we loop over i,j,k,l.

  // Row indices of non-zero entries (size = nnz):
  const int *IM = M.innerIndexPtr();
  // Values of non-zero entries (size = nnz):
  const double *valM = M.valuePtr();

  // Assuming nnz in each column is bounded by a constant, the below is O(n^2).
  // Loop over the column index j:
  for (int j = 0; j < n; j++) {
    // Loop over the non-zero entries in column j:
    for (int s = JM[j]; s < JM[j + 1]; s++) {
      int i = IM[s];         // Row index of non-zero entry #s.
      double Mij = valM[s];  // M(i,j) = Mij

      // We now compute Mij * M.
      // Loop over the column index l:
      for (int l = 0; l < n; l++) {
        // Loop over the non-zero entries in column l:
        for (int t = JM[l]; t < JM[l + 1]; t++) {
          int k = IM[t];         // Row index of non-zero entry #t.
          double Mkl = valM[t];  // M(k,l) = Mkl
          // Insert entry (k,l) of block (i,j).
          B.insert(i * n + k, j * n + l) = Mij * Mkl;
        }
      }
    }
  }
  B.makeCompressed();
  // END
 
  return B;
}
/* SAM_LISTING_END_1 */

/* @brief Solve $XA^{-1}+AX=I$ for $X$ where $A$ is s.p.d.
 * @param A, sparse s.p.d. matrix.
 * @param X, the $n\times n$ sparse matrix that solves $X+AXA=A$.
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::MatrixXd solveSpecialSylvesterEq(const Eigen::SparseMatrix<double> &A) {
  const int n = A.rows();
  assert(n == A.cols() && "Matrix A must be square.");
  Eigen::MatrixXd X(n, n);

  // TO DO: Solve the equation X*A^{-1} + A*X = I .
 
  // START
  // The equation is equivalent to X + A*X*A = A,
  // which has the system matrix C = kron(I,I) + kron(A,A).
  // Define left hand side of C*Vec(X) = b.
  Eigen::SparseMatrix<double> C(n * n, n * n);
  C.setIdentity();
  // May involve expensive copying of data?
  C += sparseKron(A);

  // Define right hand side.
  // Do not use a "sparse vector"!
  Eigen::VectorXd b(n * n);
  b.setZero();
  const int *JA = A.outerIndexPtr();
  const int *IA = A.innerIndexPtr();
  const double *valA = A.valuePtr();
  // Form b by stacking the columns of A.
  // Iterate over columns of A:
  for (int j = 0; j < n; j++) {
    // Loop over the non-zero entries in column j:
    for (int s = JA[j]; s < JA[j + 1]; s++) {
      const int i = IA[s];     // Row index of non-zero entry #s.
      b[j * n + i] = valA[s];  // b(j*n+i) = Aij
    }
  }
  // Call sparse direct solver.
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  // Since C is s.p.d., we could also use
  // SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(C);
  // x = Vec(X)
  Eigen::VectorXd x;
  x = solver.solve(b);
  X = Eigen::MatrixXd::Map(x.data(), n, n);
  // END
 
  return X;
}
/* SAM_LISTING_END_2 */
