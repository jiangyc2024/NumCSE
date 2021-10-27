#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <string>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/* SAM_LISTING_BEGIN_0 */
double extrapolate_to_pi(const unsigned int k) {
  double pi = 0.0;

  // TODO (5-8.b): Use the Aitken-Neville scheme to approximate
  // pi by extrapolation.
  // Hint: You can use the constant M_PI and Eigen's sin()
  // when calculating $s_j$ for $j=2,...,k$.

  // START

  // Data for interpolation.

  // END

  return pi;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void plotExtrapolationError(const unsigned int kmax) {
  plt::figure();
  // TODO (5-8.c): Plot the error made by extrapolate_to_pi(k) for
  // k = 2,3,...,10. Use a log-scale for the y-axis, and a linear x-axis.
  // Also, tabulate the results and errors of extrapolate_to_pi().
  // Hint: Use the constant M_PI to compute the error made by
  // extrapolate_to_pi().
  // Hint: matplotlibcpp (here = plt) implements Python's
  // matplotlib functions such as figure(), plot(), xlabel(), title(),...
  // You can use main.cpp of the LinearDataFit problem as a reference.

  // START

  // END

  std::string path = "./cx_out/pi_error_" + std::to_string(kmax) + ".png";
  plt::savefig(path);
}
/* SAM_LISTING_END_1 */
