# include <Eigen/Sparse>
# include <Eigen/SparseLU>
// or for QR: \#include <Eigen/SparseQR>
# include <Eigen/IterativeLinearSolvers> // use only if A is SPD!
/* SAM_LISTING_BEGIN_0 */
using SparseMatrix = Eigen::SparseMatrix<double>;
// Perform sparse elimination
void sparse_solve(const SparseMatrix& A, const VectorXd& b, VectorXd& x) {
  Eigen::SparseLU<SparseMatrix> solver(A);
  x = solver.solve(b);
}
/* SAM_LISTING_END_0 */
// Solve LSE, if the system matrix A is symmetric, positive definite,
// using conjugate gradient (CG) iterative solver
void sparseSpd_solve(const SparseMatrix& A, const VectorXd& b, VectorXd& x) {
  Eigen::ConjugateGradient<SparseMatrix> solver(A);
  x = solver.solve(b);
}
