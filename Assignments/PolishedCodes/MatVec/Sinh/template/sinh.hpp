#ifndef SINH_HPP
#define SINH_HPP

#include <cmath>
#include <iomanip>
#include <iostream>

/* SAM_LISTING_BEGIN_1 */
double sinh_unstable(double x) {
  const double t = std::exp(x);
  return .5 * (t - 1. / t);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void sinhError() {
  std::cout << "--> Relative error of various implementations:" << std::endl;
  std::cout << std::setw(10) << "k" << std::setw(15) << "own" << std::setw(15)
            << "std" << std::setw(15) << "err" << std::endl;

  // TODO: (1-10.a) Tabulate the relative error norms for the unstable sinh
  // function. Compare against the implementation of the standard library.
  // START

  // END
}
/* SAM_LISTING_END_2 */

#endif