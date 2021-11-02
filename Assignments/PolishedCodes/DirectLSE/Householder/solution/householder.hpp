#ifndef HOUSEHOLDER_HPP
#define HOUSEHOLDER_HPP

#include <Eigen/Dense>

/**
 * @brief
 *
 * @param x
 * @param v
 */
/* SAM_LISTING_BEGIN_0 */
void applyHouseholder(Eigen::VectorXd& x, const Eigen::VectorXd& v) {
  // TODO: (3-7.d)
  // START
  const double d = v.transpose() * x;
  x -= 2 * v * d / v.squaredNorm();
  // END
}
/* SAM_LISTING_END_0 */

/**
 * @brief
 *
 * @tparam Scalar
 * @param x
 * @param V
 */
/* SAM_LISTING_BEGIN_1 */
template <typename Scalar>
void applyHouseholder(Eigen::VectorXd& x, const Eigen::MatrixBase<Scalar>& V) {
  const unsigned int n = V.rows();

  // TODO: (3-7.e)
  // START
  for (unsigned int i = 0; i < n; ++i) {
    Eigen::VectorXd v = V.col(i);
    v(i) = std::sqrt(1. - V.col(i).squaredNorm());

    x += 2. * v.dot(x) / (1. + 2. * v.squaredNorm()) * v;
  }
  // END
}
/* SAM_LISTING_END_1 */

#endif