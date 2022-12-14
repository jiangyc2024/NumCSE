#ifndef STABRK3_H_
#define STABRK3_H_

/**
 * \file stabrk3.hpp
 * \brief NPDE homework StabRK3 code
 * \author Unknown, Oliver Rietmann, Philippe Peter
 * \date 13.04.2021
 * \copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * \brief Solves the IVP as given in task 12-7.
 *
 * \param y0 Initial state
 * \param T final time
 * \param M number of equidistant timesteps
 * \return Eigen::Vector2d state at final time
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::Vector2d PredPrey(Eigen::Vector2d y0, double T, unsigned int M) {
  const double h = T / M;
  Eigen::Vector2d y = y0;
  // TODO: (12-7.a) Solve the predator-prey model using the given RK-SSM.
  // START
  // Define right-hand-side function for Lotka-Volterra ODE
  auto f = [](Eigen::Vector2d y) -> Eigen::Vector2d {
    return {(1 - y(1)) * y(0), (y(0) - 1) * y(1)};
  };
  // Main timetepping loop: uniform stepsize
  for (unsigned int j = 0; j < M; ++j) {
    // Compute increments and updates according to \lref{def:rk} for the method
    // described by the Butcher scheme \prbeqref{eq:rkesv}
    Eigen::Vector2d k1 = f(y);
    Eigen::Vector2d k2 = f(y + h * k1);
    Eigen::Vector2d k3 = f(y + (h / 4.) * k1 + (h / 4.) * k2);
    y = y + (h / 6.) * k1 + (h / 6.) * k2 + (2. * h / 3.) * k3;
  }
  // END
  return y;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void SimulatePredPrey() {
  // TODO: (12-7.b) Simulate the predator-prey model and perform a convergence
  // study.
  // START
  // Parameters
  constexpr double T = 1.0;
  Eigen::Vector2d y0(100.0, 1.0);

  // (Approximate) reference solution
  Eigen::Vector2d y_ref = PredPrey(y0, T, std::pow(2, 14));

  Eigen::ArrayXd error(12);
  Eigen::ArrayXd M(12);
  // Studying the error for geometrically increasing numbers of equidistant
  // timesteps is the most appropriate approach to empirically exploring
  // algebraic convergence.
  M << 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192;

  // Compute errors
  for (unsigned int i = 0; i < M.size(); ++i) {
    Eigen::Vector2d y = PredPrey(y0, T, M(i));
    error(i) = (y - y_ref).norm();
  }
  // Print error table
  std::cout << std::setw(15) << "N" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;
  // Formatted output in C++
  for (unsigned int i = 0; i < M.size(); ++i) {
    std::cout << std::setw(15) << M(i) << std::setw(15) << error(i);
    if (i > 0) {
      std::cout << std::setw(15) << std::log2(error(i - 1) / error(i));
    }
    std::cout << std::endl;
  }
  // END
}
/* SAM_LISTING_END_1 */

#endif  // #ifndef STABRK3_H_
