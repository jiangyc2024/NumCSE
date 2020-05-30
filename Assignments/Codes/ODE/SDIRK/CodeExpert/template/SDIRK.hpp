#ifndef SDIRK_HPP
#define SDIRK_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "polyfit.hpp"

/* SAM_LISTING_BEGIN_0 */
//! \brief One step of autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using
//! SDIRK method Use SDIRK method for first order ode z' = f(z). Steps of size
//! h. \tparam Eigen::VectorXd type of solution space y and initial data y0
//! \param[in] z0 initial data z(0) \param[in] h size of the step \param[in]
//! gamma parameter \return next step z1
Eigen::Vector2d sdirkStep(const Eigen::Vector2d &z0, double h, double gamma) {
 
  Eigen::Vector2d res;
  // TO DO (13-3.f): compute one timestep of the ODE
  // START
  
  // END
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_3 */
//! \brief Solve autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using SDIRK
//! method Use SDIRK method for first order ode z' = f(z), with N equidistant
//! steps \tparam Eigen::VectorXd type of solution space z = [y,y']! and initial
//! data z0 = [y(0), y'(0)] \param[in] z0 initial data z(0) \param[in] N number
//! of equidistant steps \param[in] T final time of simulation \param[in] gamma
//! parameter \return vector containing each step of z_k (y and y')
std::vector<Eigen::Vector2d>
sdirkSolve(const Eigen::Vector2d &z0, unsigned int N, double T, double gamma) {
  // Solution vector
  std::vector<Eigen::Vector2d> res(N + 1);
  // TO DO (13-3.g): solve the ODE with uniform timesteps using the SDIRK method
  // START
  
  // END
  return res;
}
/* SAM_LISTING_END_3 */

double cvgSDIRK(void) {
  /* SAM_LISTING_BEGIN_9 */
  
  double conv_rate = 0;
  // TO DO (13-3.g) study the convergence rate of the method.
  // START
  
  // END
  return conv_rate;
  /* SAM_LISTING_END_9 */
}

#endif
