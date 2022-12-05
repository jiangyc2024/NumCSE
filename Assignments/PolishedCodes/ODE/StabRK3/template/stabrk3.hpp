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

  // END
  return y;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void SimulatePredPrey() {
  // TODO: (12-7.b) Simulate the predator-prey model and perform a convergence
  // study.
  // START

  // END
}
/* SAM_LISTING_END_1 */

#endif  // #ifndef STABRK3_H_
