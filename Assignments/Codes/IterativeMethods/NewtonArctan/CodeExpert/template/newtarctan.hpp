#ifndef NEWTARCTAN_HPP
#define NEWTARCTAN_HPP 

#include <cmath>
#include <iostream>
#include <limits>

/* @brief Newton's method to approximate $x^{(0)}$
 * @param[in] x0_ Initial guess
 * @param[out] x0 Final estimation of $x^{(0)}$, given convergence of Newton's method
 */
 
/* SAM_LISTING_BEGIN_0 */
double newton_arctan(double x0_ = 2) {
  // TO DO (9-3.c): define a suitable Newton iteration to find the
  // critical value x0 for F(x) = arctan(x).
  double x0 = x0_;
  // START
  
  // END
  return x0;
}
/* SAM_LISTING_END_0 */

#endif