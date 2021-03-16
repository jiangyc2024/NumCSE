#include "rkintegrator.hpp"
#include "polyfit.hpp"

//! \file rk3prey.hpp Solve prey/predator model with RK-SSM method

/* SAM_LISTING_BEGIN_0 */
double RK3prey() {
  double conv_rate = 0;
  // Dimension of state space
  unsigned int d = 2;
  // Exact value y(10) at final time T = 10 (approximated)
  Eigen::VectorXd yex(d);
  yex << 0.319465882659820, 9.730809352326228;
  
  // TO DO: Use an RKIntegrator to solve the predator/prey IVP up to T=10,
  // as described in (12-2.b). Tabulate the error as the number of steps
  // increases N = 2^7, 2^8, ..., 2^14, and estimate the convergence rate.
  // HINT: You may use polyfit() to calculate the convergence rate.
  // START
  
  
  // END
  return conv_rate;
}
/* SAM_LISTING_END_0 */
