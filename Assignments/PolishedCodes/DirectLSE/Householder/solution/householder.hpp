#ifndef HOUSEHOLDER_HPP
#define HOUSEHOLDER_HPP

#include <Eigen/Dense>

/**
 * \brief Applies a Householder transformation to x
 *
 * \param x x <- Hx
 * \param v Householder matrix
 */
/* SAM_LISTING_BEGIN_0 */
void applyHouseholder(Eigen::VectorXd& x, const Eigen::VectorXd& v) {
  // TODO: (2-7.d) Apply the Householder transformation with optimal complexity,
  // store result in x
  // START
  const double d = v.transpose() * x;
  x -= 2 * v * d / v.squaredNorm();
  // END
}
/* SAM_LISTING_END_0 */

/**
 * \brief Applies sequence of Householer transformations
 *
 * \tparam Scalar type
 * \param x input and result vector
 * \param V Householder matrices stored in columns
 */
/* SAM_LISTING_BEGIN_1 */
template <typename Scalar>
void applyHouseholder(Eigen::VectorXd& x, const Eigen::MatrixBase<Scalar>& V) {
  const unsigned int n = V.rows();

  // TODO: (2-7.e) Apply the sequence of Householder transformations to x
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