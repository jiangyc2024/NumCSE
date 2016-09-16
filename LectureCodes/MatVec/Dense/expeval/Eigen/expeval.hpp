#pragma once

#include <cmath>

using namespace std;

/* SAM_LISTING_BEGIN_0 */
double expeval(double x,
	       double tol=1e-8){
  // Initialization
  double y = 1.0, term = 1.0;
  long int k = 1;
  // \textbf{Termination} criterion
  while(abs(term) > tol*y) {
    term *= x/k;	// next summand
    y += term; // Summation
    ++k;
  }
  return y;
}
/* SAM_LISTING_END_0 */
