#pragma once

#include <cmath>

using namespace std;

/* SAM_LISTING_BEGIN_0 */
double expeval(double x, double tol=1e-8){
  // Initialization
  double y = 1., term = 1., k = 1.;
  // Termination
  while(abs(term) > tol*y){
  term *= x/k;	// next summand
  y += term; // Summation
  ++k;
  }
  return y;
}
/* SAM_LISTING_END_0 */
