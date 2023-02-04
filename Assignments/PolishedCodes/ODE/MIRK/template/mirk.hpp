#ifndef MIRK_HPP
#define MIRK_HPP

/**
 * \file mirk.hpp
 * \brief NPDE homework MIRK code
 * \author Unknown, Oliver Rietmann
 * \date 04.04.2021
 * \copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/LU>
#include <iomanip>
#include <iostream>

/**
 * \brief Perform 2 steps of the Newton method applied to F and its Jacobian DF
 *
 * \tparam Func suitable functor for F
 * \tparam Jac suitable functor for DF
 * \param F function F
 * \param DF Jacobian of function F
 * \param z initial guess
 * \return Eigen::VectorXd z after applying two steps of Newton's method
 */
/* SAM_LISTING_BEGIN_0 */
template <class Func, class Jac>
Eigen::VectorXd Newton2Steps(Func &&F, Jac &&DF, Eigen::VectorXd z) {
  // TODO: (12-6.e) Perform two steps of Newton's method.
  // START

  // END
  return z;
}
/* SAM_LISTING_END_0 */

/**
 * \brief Perform a single step of the MIRK scheme applied to the scalar ODE
 * y' = f(y)
 *
 * \tparam Func suitable functor for f
 * \tparam Jac suitable functor for Jacobian of f
 * \param f r.h.s. function
 * \param df Jacobian of f
 * \param y0 initial state
 * \param h step size
 * \return double next state
 */
/* SAM_LISTING_BEGIN_1 */
template <class Func, class Jac>
double MIRKStep(Func &&f, Jac &&df, double y0, double h) {
  // Coefficients of MIRK
  constexpr double v1 = 1.0;
  constexpr double v2 = 344.0 / 2025.0;
  constexpr double d21 = -164.0 / 2025.0;
  constexpr double b1 = 37.0 / 82.0;
  constexpr double b2 = 45.0 / 82.0;

  double y1 = y0;  // next state; overwrite this
  // TODO: (12-6.f) Perform one step of the MIRK scheme applied to the scalar
  // ODE given by f.
  // START

  // END
  return y1;
}
/* SAM_LISTING_END_1 */

/**
 * \brief Solve an ODE y' = f(y) using MIRK scheme on equidistant steps,
 * return the approximation of y(T)
 *
 * \tparam Func suitable functor for f
 * \tparam Jac suitable functor for the Jacobian of f
 * \param f r.h.s. f
 * \param df Jacobian of f
 * \param y0 initial state
 * \param T final time
 * \param M number of steps
 * \return double final state at T
 */
/* SAM_LISTING_BEGIN_2 */
template <class Func, class Jac>
double MIRKSolve(Func &&f, Jac &&df, double y0, double T, unsigned int M) {
  // Will contain next step
  double y = y0;
  // TODO: (12-6.g) Solve the ODE given by f until time T
  // START

  // END
  return y;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void cvgMIRK() {
  // TODO: (12-6.h) Tabulate the errors obtained from using the MIRK scheme on
  // the given ODE.
  // START

  // END
}
/* SAM_LISTING_END_3 */

#endif
