#pragma once

///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2019 SAM, D-MATH
/// Author(s): R. Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cassert>

/** @brief single point evaluation of Hermite polynomial
 * @param t evaluation site (between t1 and t2)
 * @param t1,t2 interval bounds
 * @param y1,y2 point values at endpoints
 * @param c1,c2 slopes at endpoints
 * @return value of cubic Hermite interpolant at t
 *
 * See Section 5.4.1 of course for details. Implements formula (5.4.1.4)
 */
inline double hermloceval(double t, double t1, double t2, double y1, double y2,
                          double c1, double c2) {
  assert((t >= t1) && (t <= t2));
  const double h = t2 - t1, a1 = y2 - y1, a2 = a1 - h * c1,
               a3 = h * c2 - a1 - a2;
  const double tau = (t - t1) / h;
  return (y1 + (a1 + (a2 + a3 * tau) * (tau - 1)) * tau);
}
