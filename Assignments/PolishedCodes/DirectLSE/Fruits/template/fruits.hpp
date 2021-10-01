#ifndef FRUITS_HPP
#define FRUITS_HPP

#include <Eigen/Dense>
#include <cassert>

/**
 * @brief Solve the LSE arising from
 * the problem description to determine the
 * prices of the fruits using Eigen
 *
 * @return Eigen::VectorXd the fruit prices
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd fruitPrice() {
  Eigen::VectorXd x = Eigen::VectorXd::Zero(6);

  // TODO: (3-1.b) Determine the fruit prices using Eigen.
  // START

  // END
  return x;
}
/* SAM_LISTING_END_1 */

#endif