#include <cmath>
#include <iostream>
#include <limits>

/* SAM_LISTING_BEGIN_1 */
double myfunction(double x) {
  double dy = std::log(2);
  double y = 0.;
  while (x > 2. * std::sqrt(2.)) {
    x /= 2.;
    y += dy;
  }  // \Label[line]{cq:1}
  while (x < std::sqrt(2.)) {
    x *= 2.;
    y -= dy;
  }  // \Label[line]{cq:2}
  double z = x - 1.;  // \Label[line]{cq:3}
  double dz = x * std::exp(-z) - 1.;
  while (std::abs(dz / z) > std::numeric_limits<double>::epsilon()) {
    z += dz;
    dz = x * std::exp(-z) - 1.;
  }
  return y + z + dz;  // \Label[line]{cq:4}
}
/* SAM_LISTING_END_1 */
