#ifndef RECURSIONORDER_HPP
#define RECURSIONORDER_HPP

#include <cmath>
#include <vector>

/**
 * @brief Guesses the maximal order of convergence, using 20 iterations in default
 * @param[in] n the number of iterations
 * @returns the order of the method
 */
/* SAM_LISTING_BEGIN_1 */
double testOrder(const unsigned int n = 20) {
  double order = 0;
  // TODO: (8-6.a) Compute the order of the method using recursion (8.6.1)
  // Use as initial error guess $e_1 = 1$ and $e_2 = 0.8$
  // START
  std::vector<double> e(n), loge(n);
  e[0] = 1; // initial error
  loge[0] = std::log(e[0]);
  e[1] = 0.8;
  loge[1] = std::log(e[1]);
  for (unsigned int k = 1; k < n - 1; ++k) {
    e[k + 1] = e[k] * std::sqrt(e[k - 1]); // error recursion
    loge[k + 1] = std::log(e[k + 1]);
  }
  order = (loge[n - 1] - loge[n - 2]) / (loge[n - 2] - loge[n - 3]);
  // END
  return order;
}
/* SAM_LISTING_END_1 */

#endif
