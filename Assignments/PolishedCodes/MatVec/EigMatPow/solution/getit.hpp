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
  // As said in the problem formulation, we may assume that the
  // following three lines have complexity $O(n^3)$.
  Eigen::EigenSolver<Eigen::MatrixXd> eig =
      Eigen::EigenSolver<Eigen::MatrixXd>(A);
  const Eigen::VectorXcd& V = eig.eigenvalues();
  const Eigen::MatrixXcd& W = eig.eigenvectors();

  // The following operation requires only a loop over
  // the dimension of $cx$, which is $n$.
  Eigen::VectorXcd cx = x.cast<std::complex<double>>();

  // The first operator* is a matrix vector multiplication
  // with complexity $O(n^2)$.
  Eigen::VectorXcd ret =
      W *
      // The componentwise power has complexity $O(n^2)$.
      // The second operator* is a vector-vector componentwise
      // multiplication, with complexity $O(n)$.
      (V.array().pow(k) *
       // In the following line, a linear system is solved,
       // operation with complexity $O(n^3)$
       (W.partialPivLu().solve(cx)).array())
          .matrix();

  return ret.real();
}
/* SAM_LISTING_END_1 */

#endif