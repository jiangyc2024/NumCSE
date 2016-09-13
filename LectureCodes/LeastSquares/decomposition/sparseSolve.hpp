# include <Eigen/Sparse>
# include <Eigen/SparseLU>
// or for QR: \#include <Eigen/SparseQR>
# include <Eigen/IterativeLinearSolvers> // use only if A is SPD!
using SparseMatrix = Eigen::SparseMatrix<double>;

void sparse_solve(const SparseMatrix& A, const VectorXd& b, VectorXd& x) {
  Eigen::SparseLU<SparseMatrix> solver(A);
  // or: Eigen::SparseQR<SparseMatrix> solver(A);
  x = solver.solve(b);
}

// solve LSE if the system matrix A is symmetric, positive definite
void sparseSpd_solve(const SparseMatrix& A, const VectorXd& b, VectorXd& x) {
  Eigen::ConjugateGradient<SparseMatrix> solver(A);
  x = solver.solve(b);
}
