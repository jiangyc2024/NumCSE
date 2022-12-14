#ifndef NLMATODE_HPP
#define NLMATODE_HPP

#include <iomanip>
#include <iostream>

#include "ode45.hpp"
#include "polyfit.hpp"

/**
 * \brief Finds an approximation of the matrix IVP $Y' = -(Y-Y')*Y$ at time $T$
 *
 * \param Y0 Initial data Y(0) (as matrix)
 * \param T final time of simulation
 * \return Eigen::MatrixXd solution at final time
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd matode(const Eigen::MatrixXd &Y0, double T) {
  Eigen::MatrixXd YT = Y0;  // overwrite this
  // TODO: (11-5.a): use the ode45 class to find an approximation
  // of the matrix IVP $Y' = -(Y-Y')*Y$ at time $T$
  // START

  // END
  return YT;
}
/* SAM_LISTING_END_1 */

/**
 * \brief Checks if invariant $Y'*Y$ is preserved
 *
 * \param M Initial data Y(0) (as matrix)
 * \param T final time of simulation
 * \return true if invariant $Y' * Y$ is preserved
 * \return false otherwise
 */
/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const Eigen::MatrixXd &M, double T) {
  bool inv_holds = false;
  // TODO: (11-5.c) check if $Y'*Y$ is preserved at the time $T$ by matode.
  // START

  // END
  return inv_holds;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
double cvgDiscreteGradientMethod() {
  double conv_rate = 0;
  // TODO: (11-5.d) compute the fitted convergence rate of the Discrete
  // gradient method. Also tabulate the values M and the errors.
  // START

  // END
  return conv_rate;
}
/* SAM_LISTING_END_3 */

#endif
