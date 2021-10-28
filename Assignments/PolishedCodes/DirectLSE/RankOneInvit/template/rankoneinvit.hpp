#ifndef RANKONEINVIT_HPP
#define RANKONEINVIT_HPP

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "timer.h"

/**
 * @brief Compute $l_{min}$ from vector $d$. Naive implementation.
 *
 * @param d An $n$-dimensional vector
 * @param tol Scalar of type 'double', the tolerance
 * @return double lmin
 */
/* SAM_LISTING_BEGIN_0 */
double rankoneinvit(const Eigen::VectorXd& d, const double& tol) {
  // Initialization
  double lmin = 0, lnew = 0;
  // TODO: (2-10.a) Port the pseudocode in algorithm (3.10.1) to C++-code.
  // START

  // END
  return lnew;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Compute $l_{min}$ from vector $d$. Optimized implementation.
 *
 * @param d An $n$-dimensional vector
 * @param tol Scalar of type 'double', the tolerance
 * @return double lmin
 */
/* SAM_LISTING_BEGIN_1 */
double rankoneinvit_fast(const Eigen::VectorXd& d, const double& tol) {
  // Initialization
  double lmin = 0, lnew = 0;
  // TODO: (2-10.c) Rewrite rankoneinvit() to have much better asymptotic
  // complexity.
  // START

  // END
  return lnew;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Tabulates the runtimes of the two different implementations for
 * different sizes of $d$.
 *
 */
/* SAM_LISTING_BEGIN_2 */
void rankoneinvit_runtime() {
  constexpr unsigned int repeats = 3;
  constexpr double tol = 1e-3;
  // TODO: (2-10.e) Tabulate the runtimes of both implementations according to
  // the problem description. The choice of the "tol" argument affects the
  // runtime.
  // START

  // END
}
/* SAM_LISTING_END_2 */

#endif
