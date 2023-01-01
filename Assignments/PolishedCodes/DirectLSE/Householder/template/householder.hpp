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

  // END
}
/* SAM_LISTING_END_1 */

#endif