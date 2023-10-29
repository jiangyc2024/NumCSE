///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>

namespace diffq {


using std::cout;
using std::setprecision;
using std::fixed;
using std::setw;
using std::abs;
using std::endl;

inline
/* SAM_LISTING_BEGIN_0 */
// Difference quotient approximation
// of the derivative of $\exp$
void diffq() {
  double h = 0.1; 
  const double x = 0.0;
  for (int i = 1; i <= 16; ++i) {
    const double df = (exp(x + h) - exp(x)) / h;
    cout << setprecision(14) << fixed;
    cout << setw(5) << -i << setw(20) << abs(df - 1) << endl;
    h /= 10;
  }
}
/* SAM_LISTING_END_0 */


} //namespace diffq