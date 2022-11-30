#ifndef ORDNOTALL_H_
#define ORDNOTALL_H_

#include <Eigen/Core>
#include <iomanip>
#include <iostream>
#include <vector>

#include "polyfit.hpp"
#include "rkintegrator.hpp"

/**
 * @brief errors Approximates the order of convergence of a scheme.
 * The scheme is a Runge Kutta scheme defined by A and b when applied
 * to the first order system y'=f(y), y(0)=y0.
 * The function ouputs error of the solutions at the time T.
 *
 * @tparam Function Type for r.h.s function f.
 * @param f The r.h.s function for the ODE.
 * @param T Final time.
 * @param y0 Initial data.
 * @param A Butcher matrix $A$.
 * @param b Butcher vector $b$.
 */
/* SAM_LISTING_BEGIN_1 */
template <class Function>
double testCvgRKSSM(const Function &f, double T, double y0,
                    const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
  // Helper object carrying out the actual explicit RK-SSM
  RKIntegrator<double> rk(A, b);
  double conv_rate = 0;
  // TODO: (11-6.a) Print the errors at each final time point and estimate the
  // convergence rate.
  // START

  // END
  return conv_rate;
}
/* SAM_LISTING_END_1 */

/**
 * @brief This function compares the convergence rates of four RK single step
 * methods: explicit Euler, trapezoidal rule, RK order 3 and classical RK
 * order 4. Comparison is done for two ODEs: 1. ODE y' = (1-y)y, y(0)=.5 and 2.
 * ODE y' = |1.1 - y| + 1, y(0)=1.
 *
 */
/* SAM_LISTING_BEGIN_2 */
void cmpCvgRKSSM() {
  // TODO: (11-6.c) Print the convergence rates of four RK single step methods
  // on the two given ODEs.
  // START

  // END
}
/* SAM_LISTING_END_2 */

#endif
