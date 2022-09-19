#ifndef SPAI_HPP
#define SPAI_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tuple>
#include <unsupported/Eigen/KroneckerProduct>
#include <vector>

using index_t = int;

/**
 * @brief Compute $B = \argmin_{X \in P(A)} |I-AX|_F$
 *
 * @param A An $n \times n$ matrix
 * @return Eigen::SparseMatrix<double> The $n \times n$ matrix $= \argmin_{X \in
 * P(A)} |I-AX|_F$
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> spai(Eigen::SparseMatrix<double>& A) {
  // Size check
  assert(A.rows() == A.cols() && "Matrix must be square!");
  unsigned int N = A.rows();

  Eigen::SparseMatrix<double> B(N, N);

  // TODO: (3-6.b) Compute the sparse approximate inverse of A by exploiting
  // Eigen's internal CCS data structure
  // START

  // END

  B.makeCompressed();
  return B;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Compute $A$
 *
 * @param n The squareroot of the matrix size wanted
 * @return Eigen::SparseMatrix<double> A $n^2 \times n^2$ sparse s.p.d. matrix
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::SparseMatrix<double> init_A(unsigned int n) {
  Eigen::SparseMatrix<double> L(n, n);
  Eigen::SparseMatrix<double> R(n, n);
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
  Eigen::SparseMatrix<double> A = kroneckerProduct(L, R);
  A.makeCompressed();
  return A;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Class implementing a symmetrized SPAI preconditioner that works with
 * Eigen's ConjugateGradient
 *
 * @tparam MatrixType a matrix
 */
/* SAM_LISTING_BEGIN_4 */
template <typename MatrixType>
class SymmetrizedSPAIPreconditioner {
 public:
  /**
   * @brief Construct a new SymmetrizedSPAIPreconditioner object
   *
   */
  SymmetrizedSPAIPreconditioner() {}

  /**
   * @brief Construct a new SymmetrizedSPAIPreconditioner object
   *
   * @param A The matrix that is the symmetrized SPAI
   */
  explicit SymmetrizedSPAIPreconditioner(const MatrixType& A) : A_(A) {}

  SymmetrizedSPAIPreconditioner& analyzePattern(const MatrixType&) {
    return *this;
  }

  SymmetrizedSPAIPreconditioner& factorize(const MatrixType&) { return *this; }

  SymmetrizedSPAIPreconditioner& compute(const MatrixType&) { return *this; }

  /**
   * @brief Applies the symmetrized SPAI to the input vector
   *
   * @tparam VectorType
   * @param b vector
   * @return const VectorType Matrix-vector product
   */
  template <typename VectorType>
  inline const VectorType solve(const VectorType& b) const {
    assert(b.size() == A_.cols());
    return A_ * b;
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }

 private:
  MatrixType A_;
};
/* SAM_LISTING_END_4 */

/**
 * @brief Compute a vector of ConjugateGradient iterations needed until
 * convergence without and with a preconditioner that is given by the SPAI
 *
 * @param L The amount of systems you want to check, system size is then
 * given by $N = (2^{1,...,L})^2$
 * @return std::vector<std::pair<unsigned int, unsigned int>> vector of tuples
 * consisting of three values: the system size, the number of iterations needed
 * without preconditioning, the number of iterations with preconditioning
 */
/* SAM_LISTING_BEGIN_3 */
std::vector<std::pair<unsigned int, unsigned int>> testSPAIPrecCG(
    unsigned int L) {
  std::vector<std::pair<unsigned int, unsigned int>> cg_iterations(L);

  // TODO: (3-6.d) Compute the iterations needed by the ConjugateGradient solver
  // without a preconditioner as well as with one

  // START

  // END

  return cg_iterations;
}
/* SAM_LISTING_END_3 */

#endif
