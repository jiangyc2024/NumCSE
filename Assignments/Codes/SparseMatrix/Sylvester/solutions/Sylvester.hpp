#include <iostream>
#include <iomanip>

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
    SparseMatrix<double> X(n,n);
    // TO DO: Fill in the entries of X.
    // Don't forget to use makeCompressed().
    // START 
    // As X is a diagonal matrix, we reserve space
    // for one non-zero entry per column.
    X.reserve( VectorXi::Constant(n,1) );
    for(int j=0; j<n; j++) {
        // X_jj = d_j / (d_j^2 + 1)
        X.insert(j,j) = diagA[j] / (diagA[j]*diagA[j] + 1.0);
    }
    X.makeCompressed();
    // END
    return X;
}
/* SAM_LISTING_END_0 */

/* @brief Compute the Kronecker product $A\otimes A$.
 * \param[in] A, an $n\times n$ sparse matrix.
 * \param[out] B, the $n^2\times n^2$ sparse matrix $B=A\otimes A$
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> sparseKron(const SparseMatrix<double> &A) {
    int n = A.rows();
    assert(n == A.cols() && "Matrix A must be square.");
    
    // A has at most maxColNNZ non-zero entries per column.
    int maxColNNZ = 5;
    
    // JA is an array of size n+1, with
    // JA[j] = number of non-zero entries in columns 0,..., (j-1) of A.
    const int *JA = A.outerIndexPtr();
    
    // Define nnzColwiseA(j) = nnz( A.col(j) )
    VectorXd nnzColwiseA(n);
    for( int j=0; j<n; j++ ) {
        nnzColwiseA(j) = JA[j+1] - JA[j];
    }
    
    // B will be the Kronecker product of A with itself.
    SparseMatrix<double> B(n*n,n*n);
    
    // TO DO: Fill in the entries of B.
    // Hint: You can use nnzColwiseA when reserving space for B.
    // Alternatively, you may assume that A has at most 5 non-zero
    // entries in each column.
    // Use A.valuePtr() and A.innerIndexPtr() to define arrays similar
    // to JA. Use those arrays to access the non-zero entries of A.
    // START
    
    VectorXd nnzColwiseB;
    
    // We have that nnz(B.col(j*n+l)) = nnz(A.col(j)) * nnz(A.col(l)).
    MatrixXd tmp = nnzColwiseA * nnzColwiseA.transpose(); // O(n^2)
    nnzColwiseB = MatrixXd::Map(tmp.data(),n*n,1);
    
    // Alternatively (assuming the columns of A have bounded nnz),
    // we can use the (less accurate) shortcut:
    // int maxCol = nnzColwiseA.maxCoeff();
    // nnzColwiseB = VectorXd::Constant(n*n,maxCol*maxCol); // O(n)
    
    // Reserve sufficient space.
    B.reserve( nnzColwiseB );
    
    // We have that B.block(i*n,j*n,n,n) = A(i,j) * A,
    // i.e. B( i*n + k, j*n + l ) = A(i,j)*A(k,l).
    // Hence, we loop over i,j,k,l.
    
    // Row indices of non-zero entries (size = nnz):
    const int *IA = A.innerIndexPtr();    
    // Values of non-zero entries (size = nnz):
    const double *valA = A.valuePtr();
    
    // Assuming nnz in each column is bounded by a constant, the below is O(n^2).
    // Loop over the column index j:
    for(int j=0; j<n; j++) {
        // Loop over the non-zero entries in column j:
        for( int s=JA[j]; s<JA[j+1]; s++ ) {
            int i = IA[s]; // Row index of non-zero entry #s.
            double Aij = valA[s]; // A(i,j) = Aij
            
            // We now compute Aij * A.
            // Loop over the column index l:
            for(int l=0; l<n; l++) {
                // Loop over the non-zero entries in column l:
                for( int t=JA[l]; t<JA[l+1]; t++) {
                    int k = IA[t]; // Row index of non-zero entry #t.
                    double Akl = valA[t]; // A(k,l) = Akl
                    // Insert entry (k,l) of block (i,j).
                    B.insert(i*n+k,j*n+l) = Aij * Akl;            
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
 * \param[in] A, sparse s.p.d. matrix.
 * \param[out] X, the $n\times n$ sparse matrix that solves $X+AXA=A$.
 */
/* SAM_LISTING_BEGIN_2 */
MatrixXd solveSpecialSylvesterEq(const SparseMatrix<double> &A) {
    int n = A.rows();
    assert(n == A.cols() && "Matrix A must be square.");
    MatrixXd X(n,n);
    
    // TO DO: Solve the equation X*A^{-1} + A*X = I .
    // START
    // The equation is equivalent to X + A*X*A = A,
    // which has the system matrix C = kron(I,I) + kron(A,A).
    // Define left hand side of C*Vec(X) = b.
    SparseMatrix<double> C(n*n,n*n);
    C.setIdentity();
    C += sparseKron( A );
    
    // Define right hand side.
    SparseMatrix<double> b(n*n,1);
    b.reserve( A.nonZeros() );
    const int *JA = A.outerIndexPtr();
    const int *IA = A.innerIndexPtr();
    const double *valA = A.valuePtr();
    // Form b by stacking the columns of A.
    // Iterate over columns of A:
    for(int j=0; j<n; j++) {
        // Loop over the non-zero entries in column j:
        for( int s=JA[j]; s<JA[j+1]; s++ ) {
            int i = IA[s]; // Row index of non-zero entry #s.
            b.insert(j*n+i,0) = valA[s]; // b(j*n+i) = Aij
        }
    }
    
    SparseLU<SparseMatrix<double>> solver;
    // Since C is s.p.d., we could also use
    // SimplicialLLT<SparseMatrix<double>> solver;
    solver.compute(C);
    // x = Vec(X)
    VectorXd x;
    x = solver.solve( b );
    X = MatrixXd::Map( x.data(), n, n );    
    // END
    
    return X;
}
/* SAM_LISTING_END_2 */
