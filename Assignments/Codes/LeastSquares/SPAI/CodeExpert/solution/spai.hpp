#ifndef SPAI_HPP
#define SPAI_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tuple>
#include <unsupported/Eigen/KroneckerProduct>
#include <vector>

using namespace Eigen;

using index_t = int;

/* @brief Compute $B = \argmin_{X \in P(A)} |I-AX|_F$
 * @param[in] A An $n \times n$ matrix
 * @param[out] B The $n \times n$ matrix $= \argmin_{X \in P(A)} |I-AX|_F$
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> spai(SparseMatrix<double>& A) {
  // Size check
  assert(A.rows() == A.cols() && "Matrix must be square!");
  unsigned int N = A.rows();

  SparseMatrix<double> B(N, N);

  // TODO: (4-6.b) Compute the sparse approximate inverse of A by exploiting
  // Eigen's internal CCS data structure
  // START

  // Needed to make sure ...Ptr functions return
  // arrays specified in CCS format
  A.makeCompressed();

  // Obtain pointers to data of A
  double* valPtr = A.valuePtr();
  index_t* innPtr = A.innerIndexPtr();
  index_t* outPtr = A.outerIndexPtr();

  // Create vector for triplets of B and reserve enough space
  std::vector<Triplet<double>> triplets;
  triplets.reserve(A.nonZeros());

  // Loop over each column of $A$
  for (unsigned int i = 0; i < N; ++i) {
    // Number of non zeros in current column of $A$
    index_t nnz_i = outPtr[i + 1] - outPtr[i];
    // If column is empty: skip column (matrix is not invertible)
    if (nnz_i == 0) continue;

    // Temporarily create a (small, dense) matrix on which normal formula
    // will be applied. We project the space $\mathcal{P}(A)$
    // onto $\mathcal{P}(a_i)$
    SparseMatrix<double> C(N, nnz_i);
    std::vector<Triplet<double>> C_triplets;
    C_triplets.reserve(nnz_i * nnz_i);

    // Need to build matrix $C$. To this end we remove all columns
    // from $A$ for which $a_i = 0$.
    // To this end: loop over all non-zero entries of the $i$-th column
    // This loop has length $n$
    for (unsigned int k = outPtr[i]; k < outPtr[i + 1]; ++k) {
      // Row of this non-zero entry
      index_t row_k = innPtr[k];
      // Number of non-zero entries for $row_k$-th column
      index_t nnz_k = outPtr[row_k + 1] - outPtr[row_k];
      // Loop over all non-zeros of $row_k$-th-column
      // This loop has length complexity $n$
      for (unsigned int l = 0; l < nnz_k; ++l) {
        unsigned int innIdx = outPtr[row_k] + l;
        C_triplets.emplace_back(
            Triplet<double>(innPtr[innIdx], k - outPtr[i], valPtr[innIdx]));
      }
    }
    C.setFromTriplets(C_triplets.begin(), C_triplets.end());
    C.makeCompressed();

    // Compute $C^\top C$ and solve normal equation.
    // Complexity of product: $O(n^3)$
    // Complexity of solve: $O(n^3)$
    // Size of $b$ is at most $n$.
    SparseMatrix<double> S = C.transpose() * C;
    MatrixXd M = MatrixXd(S);
    VectorXd xt = C.row(i).transpose();
    VectorXd b = M.partialPivLu().solve(xt);

    // Loop of length at most $n$.
    for (unsigned int k = 0; k < b.size(); ++k) {
      triplets.emplace_back(Triplet<double>(innPtr[outPtr[i] + k], i, b(k)));
    }
  }

  // Build and return SPAI preconditioner
  B.setFromTriplets(triplets.begin(), triplets.end());
  // END

  B.makeCompressed();
  return B;
}
/* SAM_LISTING_END_1 */

/* @brief Compute $A$
 * @param[in] n The squareroot of the matrix size wanted
 * @param[out] A A $n^2 \times n^2$ sparse s.p.d. matrix
 */
/* SAM_LISTING_BEGIN_2 */
SparseMatrix<double> init_A(unsigned int n) {
  SparseMatrix<double> L(n, n);
  SparseMatrix<double> R(n, n);
  L.reserve(n);
  R.reserve(n + 2 * (n - 1));
  L.insert(0, 0) = std::exp(1. / n);
  R.insert(0, 0) = 2.;
  for (unsigned int i = 1; i < n; ++i) {
    L.insert(i, i) = std::exp((i + 1.) / n);
    R.insert(i, i) = 2.;
    R.insert(i - 1, i) = -1.;
    R.insert(i, i - 1) = -1.;
  }
  SparseMatrix<double> A = kroneckerProduct(L, R);
  A.makeCompressed();
  return A;
}
/* SAM_LISTING_END_2 */

// @brief Class implementing a symmetrized SPAI preconditioner that works with
// Eigen's ConjugateGradient
/* SAM_LISTING_BEGIN_4 */
template <typename MatrixType>
class SymmetrizedSPAIPreconditioner {
 public:
  SymmetrizedSPAIPreconditioner() {}

  // @param[in] $A$ The matrix that is the symmetrized SPAI
  explicit SymmetrizedSPAIPreconditioner(const MatrixType& A) : A_(A) {}

  SymmetrizedSPAIPreconditioner& analyzePattern(const MatrixType&) {
    return *this;
  }

  SymmetrizedSPAIPreconditioner& factorize(const MatrixType&) { return *this; }

  SymmetrizedSPAIPreconditioner& compute(const MatrixType&) { return *this; }

  /* @brief Applies the symmetrized SPAI to the input vector
   * @param[in] $b$ A vector
   */
  template <typename VectorType>
  inline const VectorType solve(const VectorType& b) const {
    assert(b.size() == A_.cols());
    return A_ * b;
  }

  ComputationInfo info() { return Success; }

 private:
  MatrixType A_;
};
/* SAM_LISTING_END_4 */

/* @brief Compute a vector of ConjugateGradient iterations needed until
 * convergence without and with a preconditioner that is given by the SPAI
 * @param[in] L The amount of systems you want to check, system size is then
 * given by $N = (2^{1,...,L})^2$
 * @param[out] vector of tuples consisting of three values: the system size, the
 * number of iterations needed without preconditioning, the number of iterations
 * with preconditioning
 */
/* SAM_LISTING_BEGIN_3 */
using tuple_vector =
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>;
tuple_vector testSPAIPrecCG(unsigned int L) {
  tuple_vector cg_iterations(L);

  // TODO: (4-6.d) Compute the iterations needed by the ConjugateGradient solver
  // without a preconditioner as well as with one START initialize r.h.s. with
  // its maximal size we will use head() to get smaller versions of it
  VectorXd b = VectorXd::Ones(1 << (2 * L));

  // loop over $n = 2^{1,...,L}$, keep two running variables
  for (unsigned int n = 2, i = 0; n <= (1 << L); n <<= 1, ++i) {
    unsigned int N = n * n;
    SparseMatrix<double> A = init_A(n);
    SparseMatrix<double> B = spai(A);
    // compute Ax = b without preconditioner
    ConjugateGradient<SparseMatrix<double>, Lower | Upper,
                      IdentityPreconditioner>
        cg(A);
    VectorXd sol = cg.solve(b.head(N));

    // compute Ax = b with preconditioner
    ConjugateGradient<SparseMatrix<double>, Lower | Upper,
                      SymmetrizedSPAIPreconditioner<SparseMatrix<double>>>
        cg_preconditioner;
    // note the wrapper around B.transpose():
    // this changes the storage order which has to match when executing
    // operator+
    SymmetrizedSPAIPreconditioner<SparseMatrix<double>> preconditioner(
        0.5 * (B + SparseMatrix<double>(B.transpose())));
    cg_preconditioner.preconditioner() = preconditioner;
    cg_preconditioner.compute(A);
    VectorXd sol_prec = cg_preconditioner.solve(b.head(N));

    // check if both solutions coincide
    assert((sol - sol_prec).norm() <= 1e-9);

    cg_iterations[i] = {N, cg.iterations(), cg_preconditioner.iterations()};
  }
  // END

  return cg_iterations;
}
/* SAM_LISTING_END_3 */

#endif
