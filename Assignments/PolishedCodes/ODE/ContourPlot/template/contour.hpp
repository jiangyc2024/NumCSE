#ifndef CONTOUR_HPP
#define CONTOUR_HPP

#include <Eigen/Dense>
#include <iostream>
#include <limits>
#include <vector>

#include "ode45.hpp"

/**
 * @brief Solves an initial value problem as given in the question sheet.
 *
 * @tparam GradientFunctor Functor providing operator() returning 2d gradient
 * @param gradF the gradient
 * @param y0 initial value
 * @param T final time
 * @return Eigen::Matrix<double, 2, Eigen::Dynamic> sequence of computed states
 */
/* SAM_LISTING_BEGIN_0 */
template <typename GradientFunctor>
Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePoints(
    GradientFunctor &&gradF, Eigen::Vector2d y0, double T) {
  Eigen::Matrix<double, 2, Eigen::Dynamic> States =
      Eigen::Matrix2Xd::Zero(2, 1);
  // TODO: (11-9.b) Solve the IVP as given in the task description.
  // START

  // END
  return States;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Computes the crooked egg curve.
 *
 * @return Eigen::Matrix<double, 2, Eigen::Dynamic> discretized curve
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEgg() {
  // crookedEggCurve will need to be reshaped to 2*(N+1).
  Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEggCurve =
      Eigen::Matrix2Xd::Zero(2, 1);
  // TODO: (11-9.c) Compute the crooked egg curve using computeIsolinePoints().
  // START

  // END
  return crookedEggCurve;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solves an initial value problem as given in the question sheet using
 * the gradient obtained by symmetric difference quotients.
 *
 * @tparam FFunctor functor providing operator() returning the value of F
 * @param F the function
 * @param y0 initial value
 * @param T final time
 * @return Eigen::Matrix<double, 2, Eigen::Dynamic> sequence of computed states
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FFunctor>
Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePointsDQ(
    FFunctor &&F, Eigen::Vector2d y0, double T) {
  // States will need to be reshaped to 2*(N+1).
  Eigen::Matrix<double, 2, Eigen::Dynamic> States =
      Eigen::Matrix2Xd::Zero(2, 1);
  // TODO: (11-9.d) Solve the IVP as given in the task description by using the
  // symmtric difference quotient for the gradient of F. You may copy code from
  // above.
  // START

  // END
  return States;
}
/* SAM_LISTING_END_2 */

#endif