#ifndef SDIRK_HPP
#define SDIRK_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "polyfit.hpp"

/**
 * \brief One step of autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using
 * SDIRK method Use SDIRK method for first order ode z' = f(z). Steps of size
 * h.
 *
 * \param z0 initial data z(0)
 * \param h size of the step
 * \param gamma parameter
 * \return Eigen::Vector2d next step z1
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::Vector2d sdirkStep(const Eigen::Vector2d &z0, double h, double gamma) {
  Eigen::Vector2d res = z0;
  // TODO: (12-3.f) compute one timestep of the ODE
  // START

  // END
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_9 */
double cvgSDIRK() {
  double conv_rate = 0.;
  // Initial data z0 = [y(0), y'(0)]
  Eigen::Vector2d z0;
  z0 << 1, 0;
  // Final time
  constexpr double T = 10;
  // Parameter
  const double gamma = (3. + std::sqrt(3.)) / 6.;
  // Mesh sizes
  Eigen::ArrayXd err(10);
  Eigen::ArrayXd N(10);
  N << 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240;
  // TODO: (12-3.g) Study the convergence rate of the method.
  // START

  // END
  return conv_rate;
}
/* SAM_LISTING_END_9 */

#endif
