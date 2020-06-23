#ifndef QUADINF_HPP
#define QUADINF_HPP

#include <iostream>
#include <cmath>

#include "golubwelsh.hpp"
#include "matplotlibcpp.h"

#define PI M_PI

namespace plt = matplotlibcpp;

// TO DO (8-3.e): define the function quadinf that integrates the input 
// lambda f over the real axis. First trasform the integrand with a change
// of variables and then use an n point Gauss quadrature.
// Hint 1: lambda functions can take parameters inside the [] brackets
// Hint 2: you may write an auxiliary function to compute the quadrature over 
//         a bounded interval.

//START

/* SAM_LISTING_BEGIN_1 */

// Optional: define an auxiliary function here as follows
/*
template <class Function> double quad(...) {

}
*/

/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <class Function>
double quadinf(const int n, Function &&f) {
  return NAN;
}
/* SAM_LISTING_END_2 */

//END



//! @brief perform convergence test for h(t) := exp(-(t-1)^2) and plot error

/* SAM_LISTING_BEGIN_3 */
void cvgQuadInf(void) {
    // Number of max Gauss pts.
    const int N = 100;
    plt::figure();
    //TO DO (8-3.f): plot convergence errors against number of integration nodes
    //START
    
    //END
    plt::savefig("./cx_out/convergence.png");
}
/* SAM_LISTING_BEGIN_3 */

#endif
