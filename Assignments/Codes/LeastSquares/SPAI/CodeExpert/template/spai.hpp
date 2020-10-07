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
    std::vector<std::tuple<unsigned int, unsigned int>>;
tuple_vector testSPAIPrecCG(unsigned int L) {
  tuple_vector cg_iterations(L);

  // TODO: (4-6.d) Compute the iterations needed by the ConjugateGradient solver
  // without a preconditioner as well as with one 
  // START

  // END

  return cg_iterations;
}
/* SAM_LISTING_END_3 */

#endif
