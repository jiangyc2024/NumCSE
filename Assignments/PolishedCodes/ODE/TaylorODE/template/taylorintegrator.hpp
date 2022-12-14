#ifndef TAYLORINTEGRATORHPP
#define TAYLORINTEGRATORHPP

#include <Eigen/Dense>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "ode45.hpp"

/**
 * \brief Solves the Predator Prey model using the Taylor expansion method.
 *
 * \param alpha1 parameter of ODE
 * \param beta1 parameter of ODE
 * \param alpha2 parameter of ODE
 * \param beta2 parameter of ODE
 * \param T final time
 * \param y0 initial state
 * \param M number of steps to take
 * \return std::vector<Eigen::Vector2d> of states
 */
/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Vector2d> SolvePredPreyTaylor(double alpha1, double beta1,
                                                 double alpha2, double beta2,
                                                 double T,
                                                 const Eigen::Vector2d &y0,
                                                 unsigned int M) {
  // Vector storing the states
  std::vector<Eigen::Vector2d> res;
  // TODO: (11-8.c) Solve the predator prey model using the Taylor expansion
  // method.
  // START

  // END
  return res;
}
/* SAM_LISTING_END_1 */

void PrintErrorTable(const Eigen::ArrayXd &M, const Eigen::ArrayXd &error) {
  std::cout << std::setw(15) << "M" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;

  for (unsigned int i = 0; i < M.size(); ++i) {
    std::cout << std::setw(15) << M(i) << std::setw(15) << error(i);
    if (i > 0) {
      std::cout << std::setw(15) << std::log2(error(i - 1) / error(i));
    }
    std::cout << std::endl;
  }
}

/* SAM_LISTING_BEGIN_2 */
double testCvgTaylorMethod() {
  double cvgrate = 0;
  // TODO: (11-8.d) Generate an error table and return an estimate of the
  // convergence rate of the Taylor method.
  // START

  // END
  return cvgrate;
}
/* SAM_LISTING_END_2 */

#endif
