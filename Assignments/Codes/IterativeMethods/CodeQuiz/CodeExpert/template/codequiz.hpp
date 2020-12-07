#pragma once

#include <cmath>
#include <iostream>
#include <limits>

/* SAM_LISTING_BEGIN_1 */
double myfunction(double x) {
  double log2 = 0.693147180559945;
  double y = 0.;
  while (x > 2. * std::sqrt(2.)) {
    x /= 2.;
    y += log2;
  }  // \Label[line]{cq:1}
  while (x < std::sqrt(2.)) {
    x *= 2.;
    y -= log2;
  }                   // \Label[line]{cq:2}
  double z = x - 1.;  // \Label[line]{cq:3}
  double dz = x * std::exp(-z) - 1.;
  while (std::abs(dz / z) > std::numeric_limits<double>::epsilon()) {
    z += dz;
    dz = x * std::exp(-z) - 1.;
  }
  return y + z + dz;  // \Label[line]{cq:4}
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double myfunction_modified(double x) {
  double log2 = 0.693147180559945;
  double y = 0.;
  while (x > 2. * std::sqrt(2.)) {
    x /= 2.;
    y += log2;
  }
  while (x < std::sqrt(2.)) {
    x *= 2.;
    y -= log2;
  }
  double z = x - 1.;
  double dz = x * std::exp(-z) - 1.;
  // TO DO: (9-2.e) Write a for-loop that achieves the same accuracy
  // as the third while-loop of myfunction(). I.e. fix the number of
  // iterations before looping.
  // START

  // END
  return y + z + dz;
}
/* SAM_LISTING_END_2 */
