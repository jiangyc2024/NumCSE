#ifndef DISCRETECONVOLUTION_HPP
#define DISCRETECONVOLUTION_HPP

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <unsupported/Eigen/FFT>

/**
 * @brief A ``naive'' two-loop implementation of discrete periodic convolution
 *
 * @param u vector of size n
 * @param x vector of size n to be convolved with u
 * @return Eigen::VectorXd discrete periodic convolution of u and x
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd pconv(const Eigen::VectorXd &u, const Eigen::VectorXd &x) {
  const int n = x.size();
  assert(n == u.size());
  Eigen::VectorXd z = Eigen::VectorXd::Zero(n);
  for (int k = 0; k < n; ++k) {
    // Calculate z[k] = \sum_{j=0}^{n-1} u[(k-j)%n] * x[j]
    for (int j = 0, l = k; j <= k; ++j, --l) {
      z[k] += u[l] * x[j];
    }
    for (int j = k + 1, l = n - 1; j < n; ++j, --l) {
      z[k] += u[l] * x[j];
    }
  }
  return z;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Efficient discrete periodic convolution
 *
 * @param p vector of size n
 * @param x vector of size n to be convolved with p
 * @return Eigen::VectorXd discrete periodic convolution of p and x
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd pconv_fast(const Eigen::VectorXd &p, const Eigen::VectorXd &x) {
  const int n = x.size();
  assert(n == p.size());
  Eigen::VectorXd z = Eigen::VectorXd::Zero(n);

  // TODO: (4-5.b) Compute the discrete periodic convolution of p with x
  // in an efficient manner.
  // START

  // END
  return z;
}
/* SAM_LISTING_END_2 */

#endif