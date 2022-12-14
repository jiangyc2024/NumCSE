#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>

#include "polyfit.hpp"

/**
 * @brief Performs a single step of classical Runge-Kutta method of order 4.
 *
 * @tparam Function providing operator()
 * @tparam State vector overloaded with basic operations
 * @param odefun The function describing the ODE
 * @param h step size
 * @param y0 initial state
 * @param y1 next state
 */
/* SAM_LISTING_BEGIN_1 */
template <class Function, class State>
void rk4step(Function &&odefun, double h, const State &y0, State &y1) {
  // TODO: (11-4.c) Implement a single step of RK4 for
  // the ODE y' = odefun(y), starting from y0 and using
  // the step size h. Save the result in y1.
  // START

  // END
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_0 */
double testcvgRK4() {
  // Parameters and initial values:
  constexpr double T = 1;
  constexpr unsigned int n = 5;
  Eigen::VectorXd y0(2 * n);
  for (unsigned int i = 0; i < n; ++i) {
    y0(i) = (i + 1.) / n;
    y0(i + n) = -1;
  }

  // TODO: (11-4.d) Implement the function f from (12-4.a) as a lambda
  // function, i.e. the right hand side of the system of first order ODEs.
  // START

  // END

  double conv_rate = 0;
  // Table header
  std::cout << std::setw(8) << "N" << std::setw(20) << "Error" << std::endl;

  // TODO: (11-4.d) Tabulate the error at time T=1, using
  // N=2,4,...,1024 RK4 steps. Use $N=2^{12}$ steps to calculate
  // the "exact" solution. Then, estimate the algebraic convergence
  // rate of the errors.
  // HINT: You can use polyfit() to calculate the convergence rate.
  // START

  // END

  return conv_rate;
}
/* SAM_LISTING_END_0 */

#endif