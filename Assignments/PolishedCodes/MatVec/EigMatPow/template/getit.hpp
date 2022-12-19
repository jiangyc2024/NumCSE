#ifndef GETIT_HPP
#define GETIT_HPP

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

/**
 * \brief Performs the computation $y = A^k x$.
 *
 * \param A matrix
 * \param x vector
 * \param k exponent
 * \return Eigen::VectorXd $y = A^k x$
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd getit(const Eigen::MatrixXd& A, const Eigen::VectorXd& x,
                      unsigned int k) {
  Eigen::EigenSolver<Eigen::MatrixXd> eig =
      Eigen::EigenSolver<Eigen::MatrixXd>(A);
  const Eigen::VectorXcd& V = eig.eigenvalues();
  const Eigen::MatrixXcd& W = eig.eigenvectors();

  Eigen::VectorXcd cx = x.cast<std::complex<double>>();

  Eigen::VectorXcd ret =
      W * (V.array().pow(k) * (W.partialPivLu().solve(cx)).array()).matrix();

  return ret.real();
}
/* SAM_LISTING_END_1 */

#endif