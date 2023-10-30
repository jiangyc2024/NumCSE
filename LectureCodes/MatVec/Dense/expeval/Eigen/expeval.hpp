///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>
#include <cstdint>

namespace expeval {


using std::abs;

inline
/* SAM_LISTING_BEGIN_0 */
double expeval(double x,
	       double tol=1e-8){
  // Initialization
  double y = 1.0;
  double term = 1.0;
  uint64_t k = 1;
  // \textbf{Termination} criterion
  while(abs(term) > tol*y) {
    term *= x/static_cast<double>(k);	// next summand
    y += term; // Summation
    ++k;
  }
  return y;
}
/* SAM_LISTING_END_0 */


} //namespace expeval