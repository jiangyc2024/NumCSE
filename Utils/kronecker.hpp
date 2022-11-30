#ifndef NCSE_KRONECKER_HPP
#define NCSE_KRONECKER_HPP

#include <Eigen/Dense>

inline Eigen::MatrixXd kron(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
  const auto dimensions = std::make_tuple(A.rows(), A.cols(), B.rows(), B.cols());
  const auto [rA, cA, rB, cB] = dimensions;
  Eigen::MatrixXd K(rA*rB, cA*cB);
  for (Eigen::Index i = 0; i < rA; ++i) {
    for (Eigen::Index j = 0; j < cA; ++j) {
      K.block(i*rB, j*cB, rB, cB) = A(i,j)*B;
    }
  }
  return K;
}

# endif // NCSE_KRONECKER_HPP
