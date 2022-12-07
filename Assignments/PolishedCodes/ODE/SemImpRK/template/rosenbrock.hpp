#ifndef ROSENBROCK_HPP
#define ROSENBROCK_HPP

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#include "polyfit.hpp"

/**
 * \brief Solve the autonomous IVP y' = f(y), y(0) = y0 using Rosenbrock method
 * Use semi-implicit Rosenbrock method using Jacobian evaluation. Equidistant
 * steps of size T/N.
 *
 * \tparam Function function type for r.h.s. f
 * \tparam Jacobian function type for Jacobian df
 * \param f r.h.s. func f
 * \param df Jacobian df of f
 * \param y0 initial data y(0)
 * \param N number of equidistant steps
 * \param T final time
 * \return std::vector<Eigen::VectorXd> of y_k for each step k from 0 to N
 */
/* SAM_LISTING_BEGIN_0 */
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> solveRosenbrock(Function &&f, Jacobian &&df,
                                             const Eigen::VectorXd &y0,
                                             unsigned int N, double T) {
  // Will contain all time steps
  std::vector<Eigen::VectorXd> res(N + 1);

  // TODO: (12-4.e) Implement the Rosenbrock method. Note that the function
  // f should take as argument an Eigen::VectorXd, and return one as well.
  // The Jacobian df should take as argument an Eigen::VectorXd, but return
  // a square Eigen::MatrixXd.
  // START

  // END
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
double cvgRosenbrock() {
  double cvgRate = 0;
  // TODO: (12-4.f) Use polyfit() to estimate the rate of convergence
  // for solveRosenbrock().
  // START

  // END
  return cvgRate;
}
/* SAM_LISTING_END_1 */

#endif