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
  // Right-hand-side vector field of isoline ODE
  auto rhs = [gradF](const Eigen::Vector2d &x) {
    Eigen::Vector2d gradFx = gradF(x);
    Eigen::Vector2d out(-gradFx(1), gradFx(0));  // rotate
    out /= gradFx.norm();                        // scale
    return out;
  };
  // Adaptive explicit Runge-Kutta method
  ode45<Eigen::Vector2d> integrator(rhs);
  // Set tolerances for timestep control
  integrator.options.rtol = 0.00001;
  integrator.options.atol = 1e-9;
  // Perform explicit timestepping
  std::vector<std::pair<Eigen::Vector2d, double>> sol = integrator.solve(y0, T);
  // Convert output into requested format: points on isoline arranged into the
  // columns of a matrix.
  const unsigned int N = sol.size() - 1;
  States = Eigen::MatrixXd::Zero(2, N + 1);
  for (unsigned int i = 0; i <= N; i++) {
    States.col(i) = sol[i].first;
  }
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
  Eigen::Vector2d y0(1., 0.);
  constexpr double T = 4.;
  auto gradF = [](const Eigen::Vector2d &x) {
    double x2 = x(0) * x(0), y2 = x(1) * x(1);
    double x2y2 = 4. * (x2 + y2);
    Eigen::Vector2d grad(x2y2 * x(0) - 3. * x2, x2y2 * x(1) - 3. * y2);
    return grad;
  };
  crookedEggCurve = computeIsolinePoints(gradF, y0, T);
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
  // Get gradient of F by symmetric difference quotients with a
  // span given by the root of the machine precision in order to curb
  // cancellation as much as possible.
  const double epsrt = std::sqrt(std::numeric_limits<double>::epsilon());
  auto gradF = [epsrt, F](const Eigen::Vector2d &x) {
    Eigen::Vector2d grad, dx, dy;
    dx << epsrt, 0.;
    dy << 0., epsrt;
    grad(0) = (F(x + dx) - F(x - dx)) / (2 * epsrt);
    grad(1) = (F(x + dy) - F(x - dy)) / (2 * epsrt);
    return grad;
  };
  // The remainder as in computeIsolinePoints()
  auto rhs = [gradF](const Eigen::Vector2d &x) {
    Eigen::Vector2d gradFx = gradF(x);
    Eigen::Vector2d out(-gradFx(1), gradFx(0));  // rotate
    out /= gradFx.norm();                        // scale
    return out;
  };

  ode45<Eigen::Vector2d> integrator(rhs);
  integrator.options.rtol = 0.00001;
  integrator.options.atol = 1e-9;

  std::vector<std::pair<Eigen::Vector2d, double>> sol = integrator.solve(y0, T);
  const unsigned int N = sol.size() - 1;
  States = Eigen::MatrixXd::Zero(2, N + 1);
  for (unsigned int i = 0; i <= N; i++) {
    States.col(i) = sol[i].first;
  }
  // END
  return States;
}
/* SAM_LISTING_END_2 */

#endif