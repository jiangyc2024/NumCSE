#ifndef NEWTARCTAN_HPP
#define NEWTARCTAN_HPP

#include <cmath>
#include <iostream>
#include <limits>

/* @brief Newton's method to approximate $x^{(0)}$
 * @param[in] x0_ Initial guess
 * @param[out] x0 Final estimation of $x^{(0)}$, given convergence of Newton's
 * method
 */

/* SAM_LISTING_BEGIN_0 */
double newton_arctan(double x0_ = 2) {
  // TO DO (9-3.c): define a suitable Newton iteration to find the
  // critical value x0 for F(x) = arctan(x).
  double x0 = x0_;
  // START
  double upd = 1;
  double eps = std::numeric_limits<double>::epsilon();

  while (upd > eps) {
    // Newton iteration for g(x)
    double x1 =
        (-x0 + (1 - x0 * x0) * std::atan(x0)) / (1 - 2 * x0 * std::atan(x0));
    // Relative size of the Newton update
    upd = std::abs((x1 - x0) / x1);
    x0 = x1;
  }
  // END
  return x0;
}
/* SAM_LISTING_END_0 */

#endif
