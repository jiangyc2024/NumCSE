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
  // Matrix A for evaluation of f
  Eigen::Matrix2d A;
  A << 0., 1., -1., -1.;
  // Reuse factorization
  auto A_lu = (Eigen::Matrix2d::Identity() - h * gamma * A).partialPivLu();
  Eigen::Vector2d az = A * z0;

  // Increments
  Eigen::Vector2d k1 = A_lu.solve(az);
  Eigen::Vector2d k2 = A_lu.solve(az + h * (1 - 2 * gamma) * A * k1);

  // Next step
  res = z0 + h * 0.5 * (k1 + k2);
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

  // Exact solution (only y(t)) given z0 = [y(0), y'(0)] and t
  auto yex = [&z0](double t) {
    return 1. / 3. * std::exp(-t / 2.) *
           (3. * z0(0) * std::cos(std::sqrt(3.) * t / 2.) +
            std::sqrt(3.) * z0(0) * std::sin(std::sqrt(3.) * t / 2.) +
            2. * std::sqrt(3.) * z0(1) * std::sin(std::sqrt(3.) * t / 2.));
  };

  // Store old error for rate computation
  double errold = 0;
  std::cout << std::setw(15) << "n" << std::setw(15) << "maxerr"
            << std::setw(15) << "rate" << std::endl;
  // Loop over all meshes
  for (unsigned int i = 0; i < N.size(); ++i) {
    unsigned int n = N(i);
    // Get solution
    std::vector<Eigen::Vector2d> sol(n + 1);
    // Equidistant step size
    const double h = T / n;
    // Push initial data
    sol.at(0) = z0;
    // Main loop
    for (unsigned int i = 1; i <= n; ++i) {
      sol.at(i) = sdirkStep(sol.at(i - 1), h, gamma);
    }
    // Compute error
    err(i) = std::abs(sol.back()(0) - yex(T));

    // Print table
    std::cout << std::setw(15) << n << std::setw(15) << err(i);
    if (i > 0) std::cout << std::setw(15) << std::log2(errold / err(i));
    std::cout << std::endl;

    // Store old error
    errold = err(i);
  }
  Eigen::VectorXd coeffs = polyfit(N.log(), err.log(), 1);
  conv_rate = -coeffs(0);
  // END
  return conv_rate;
}
/* SAM_LISTING_END_9 */

#endif
