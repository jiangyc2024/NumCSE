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
  
  // END
}
/* SAM_LISTING_END_1 */

#endif
