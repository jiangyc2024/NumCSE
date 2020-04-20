// TO DO: (3-10.a) Write the function rankoneinvit().
// (3-10.c) Write the function rankoneinvit_fast().
// (3-10.e) Write a function "void rankoneinvit_runtime()" that tabulates
// the runtimes of both implementations according to the problem description.
// This subtask can be solved in many ways, and is not a part of the tests.
// In particular, the choice of the "tol" argument affects the runtime.
#ifndef RANKONEINVIT_HPP
#define RANKONEINVIT_HPP

#include "timer.h"
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

using namespace Eigen;

/* @brief Compute $l_{min}$ from vector $d$
 * Naive implementation
 * @param[in] d An $n$-dimensional vector
 * @param[in] tol Scalar of type 'double', the tolerance
 * @param[out] lmin Scalar of type 'double'
 */
/* SAM_LISTING_BEGIN_0 */
double rankoneinvit(const VectorXd &d, const double &tol) {
  // START
  
  // placeholder for compilation:
  return 0.;
  
  // END
}
/* SAM_LISTING_END_0 */

/* @brief Compute $l_{min}$ from vector $d$
 * Optimized implementation
 * @param[in] d An $n$-dimensional vector
 * @param[in] tol Scalar of type 'double', the tolerance
 * @param[out] lmin Scalar of type 'double'
 */
/* SAM_LISTING_BEGIN_1 */
double rankoneinvit_fast(const VectorXd &d, const double &tol) {
  // START
  
  // placeholder for compilation:
  return 0.;
    
  // END
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void rankoneinvit_runtime() {
  // START
  
  // END
}
/* SAM_LISTING_END_2 */

#endif
