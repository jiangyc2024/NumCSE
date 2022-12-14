#ifndef RK3PREY_HPP
#define RK3PREY_HPP

#include "polyfit.hpp"
#include "rkintegrator.hpp"

/* SAM_LISTING_BEGIN_0 */
double RK3prey() {
  double conv_rate = 0;
  // Dimension of state space
  constexpr unsigned int d = 2;
  // Exact value y(10) at final time T = 10 (approximated)
  Eigen::VectorXd yex(d);
  yex << 0.319465882659820, 9.730809352326228;

  // Implementation of butcher scheme
  constexpr unsigned int s = 3;
  Eigen::MatrixXd A(s, s);
  Eigen::VectorXd b(s);
  A << 0, 0, 0, 1. / 3., 0, 0, 0, 2. / 3., 0;
  b << 1. / 4., 0, 3. / 4.;

  // TODO: (11-2.b) Use an RKIntegrator to solve the predator/prey IVP up to
  // T=10, as described in the task description. Tabulate the error as the
  // number of steps increases N = 2^7, 2^8, ..., 2^14, and estimate the
  // convergence rate.
  // HINT: You may use polyfit() to calculate the convergence
  // rate.
  // START

  // END
  return conv_rate;
}
/* SAM_LISTING_END_0 */

#endif