#ifndef NLMATODE_HPP
#define NLMATODE_HPP

#include "ode45.hpp"
#include "polyfit.hpp"
using namespace Eigen;


//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
/* SAM_LISTING_BEGIN_1 */
MatrixXd matode(const MatrixXd &Y0, double T) {
  // TO DO (12-5.a): use the ode45 class to find an approximation
  // of the matrix IVP $Y' = -(Y-Y')*Y$ at time $T$
  MatrixXd YT;
  // START
  
  // END
  return YT;
}
/* SAM_LISTING_END_1 */


//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const MatrixXd &M, double T) {
  // TO DO (12-5.c): check if $Y'*Y$ is preserved at the time $T$ by matode.
  // START
  return false;
  // END
}
/* SAM_LISTING_END_2 */


/* SAM_LISTING_BEGIN_3 */
double cvgDiscreteGradientMethod(void) {
  // TO DO (12-5.d): compute the fitted convergence rate of the Discrete gradient
  // method. Also tabulate the values M and the errors. 
  double conv_rate;
  // START
  
  // END
  return conv_rate;
}
/* SAM_LISTING_END_3 */

#endif