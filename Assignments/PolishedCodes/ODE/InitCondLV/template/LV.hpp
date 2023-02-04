#ifndef LV_HPP
#define LV_HPP
////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>

#include "ode45.hpp"

/**
 * @brief Compute the maps Phi and W at time T.
 * Use initial data given by u0 and v0.
 *
 * @param u0 First component.
 * @param v0 Second component.
 * @param T Final time.
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<Eigen::Vector2d, Eigen::Matrix2d> PhiAndW(double u0, double v0,
                                                    double T) {
  std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW;
  // TODO: (11-7.f)  Calculate the time evolution of Phi(t,y0) and W(t,y0)
  // up to the final time T, using the inital value y0=[u0,v0].
  // Save the values of Phi and W at time T in PaW.first and
  // PaW.second respectively.
  // Hint: First write the ODE
  // [ y'(t), W'(t,y0) ] = [ f(y(t)), Df(y(t))W(t,y0) ]
  // in vectorized form, i.e. find a function g:R^6 -> R^6
  // such that the ODE w'=g(w) is equivalent.
  // Then use the ode45 class for the function g.
  // START

  // END
  return PaW;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Compute initial conditions u0, v0 such that solution has period T=5
 * inital approximation [3,2]^T
 *
 * @return std::pair<double, double> initial conditions
 */
/* SAM_LISTING_BEGIN_2 */
std::pair<double, double> findInitCond() {
  Eigen::Vector2d y;

  // TODO: (11-7.g) Use Newton's method with initial guess (3,2) to find a
  // zero of the function F from (12-7.c) for the period $T_P=5$.
  // Store the resulting vector in y.
  // Hint: F and its Jacobian matrix can be computed from the result of
  // PhiAndW().
  // START

  // END

  // TODO: (11-7.g) Calculate the evolution of the Lotka-Volterra ODE up to
  // time 100, using the initial value y, found by the Newton iterations above.
  // If the evolution is indeed 5-periodic, then the solution at time 100
  // should have the value y. Print a warning message if this is not the case.
  // START

  // END

  return std::make_pair(y(0), y(1));
}
/* SAM_LISTING_END_2 */

#endif
