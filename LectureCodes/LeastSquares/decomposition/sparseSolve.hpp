#include <Eigen/Sparse>

#include <Eigen/SparseLU>

// or for QR: \#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers> // use only if A is SPD!

#include <stdexcept>

namespace sparseSolve {


using Eigen::MatrixXd;
using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
using SparseMatrix = Eigen::SparseMatrix<double>;
// Perform sparse elimination
inline void sparse_solve(const SparseMatrix &A, const VectorXd &b, VectorXd &x) {
  const Eigen::SparseLU<SparseMatrix> solver(A);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Matrix factorization failed");
  }
  x = solver.solve(b);
}
/* SAM_LISTING_END_0 */

// Solve LSE, if the system matrix A is symmetric, positive definite,
// using conjugate gradient (CG) iterative solver
inline void sparseSpd_solve(const SparseMatrix &A, const VectorXd &b, VectorXd &x) {
  const Eigen::ConjugateGradient<SparseMatrix> solver(A);
  x = solver.solve(b);
}


} //namespace sparseSolve