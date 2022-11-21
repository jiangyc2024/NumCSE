#ifndef SDIRK_HPP
#define SDIRK_HPP

#include "polyfit.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

//! \brief One step of autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using
//! SDIRK method Use SDIRK method for first order ode z' = f(z). Steps of size
//! h.
//! \tparam Eigen::VectorXd type of solution space y and initial data y0
//! \param[in] z0 initial data z(0)
//! \param[in] h size of the step
//! \param[in] gamma parameter
//! \return next step z1
/* SAM_LISTING_BEGIN_0 */
Eigen::Vector2d sdirkStep(const Eigen::Vector2d &z0, double h, double gamma) {

  Eigen::Vector2d res;
  // TO DO (12-3.f): compute one timestep of the ODE
  // START

  // END
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_9 */
double cvgSDIRK() {

  double conv_rate;
  // TO DO (12-3.g) study the convergence rate of the method.
  // START

  // END
  return conv_rate;
  /* SAM_LISTING_END_9 */
}

#endif
