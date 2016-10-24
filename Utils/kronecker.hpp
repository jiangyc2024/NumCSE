# ifndef NCSE_KRONECKER_HPP
# define NCSE_KRONECKER_HPP

# include <Eigen/Dense>
using Eigen::MatrixXd;

MatrixXd kron(const MatrixXd& A, const MatrixXd& B) {
  const long rA = A.rows(), cA = A.cols(),
             rB = B.rows(), cB = B.cols();
  MatrixXd K(rA*rB, cA*cB);
  for (long i = 0; i < rA; ++i) {
    for (long j = 0; j < cA; ++j) {
      K.block(i*rB, j*cB, rB, cB) = A(i,j)*B;
    }
  }
  return K;
}

# endif // NCSE_KRONECKER_HPP
