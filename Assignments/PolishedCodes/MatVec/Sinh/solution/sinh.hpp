#ifndef SINH_HPP
#define SINH_HPP

#include <cmath>
#include <iomanip>
#include <iostream>

#include "sinh_stable.hpp"

/* SAM_LISTING_BEGIN_1 */
double sinh_unstable(double x) {
  double t = std::exp(x);
  return .5 * (t - 1. / t);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void sinhError() {}
/* SAM_LISTING_END_2 */

#endif