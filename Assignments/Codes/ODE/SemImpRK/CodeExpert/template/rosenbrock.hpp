#include <Eigen/Dense>

#include <iostream>
#include <iomanip>
#include <vector>
#include "polyfit.hpp"

using namespace Eigen;

//! \brief Solve the autonomous IVP y' = f(y), y(0) = y0 using Rosenbrock method
//! Use semi-implicit Rosenbrock method using Jacobian evaluation. Equidistant steps of size T/N.
//! \tparam Function function type for r.h.s. f
//! \tparam Jacobian function type for Jacobian df
//! \param[in] f r.h.s. func f
//! \param[in] df Jacobian df of f
//! \param[in] y0 initial data y(0)
//! \param[in] N number of equidistant steps
//! \param[in] T final time
//! \return vector of y_k for each step k from 0 to N
/* SAM_LISTING_BEGIN_0 */
template <class Function, class Jacobian>
std::vector<VectorXd> solveRosenbrock(Function &&f, Jacobian &&df,
                                       const VectorXd &y0, unsigned int N, double T) {
  
  // Will contain all time steps
  std::vector<VectorXd> res(N+1);
  
  // TO DO: (13-2.c) Implement the Rosenbrock method. Note that the function
  // f should take as argument an Eigen::VectorXd, and return one as well.
  // The Jacobian df should take as argument an Eigen::VectorXd, but return
  // a square Eigen::MatrixXd.
  // START
  
  // END
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
double cvgRosenbrock(void) {
  
  double cvgRate=0;
  // TO DO: (13-2.d) Use polyfit() to estimate the rate of convergence
  // for solveRosenbrock().
  // START
  
  // END
  return cvgRate;
}
/* SAM_LISTING_END_1 */
