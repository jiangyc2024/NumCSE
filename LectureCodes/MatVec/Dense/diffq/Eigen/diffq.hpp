#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

/* SAM_LISTING_BEGIN_0 */
//! Difference quotient approximation
//! of the derivative of $\exp$
void diffq(){
  double h = 0.1, x = 0.0;
  for(int i = 1; i <= 16; ++i){
    double df = (exp(x+h)-exp(x))/h;
    cout << setprecision(14) << fixed;
    cout  << setw(5) <<  -i
	  << setw(20) << df-1 << endl;
    h /= 10;
  }
}
/* SAM_LISTING_END_0 */
