#ifndef CROSS_HPP
#define CROSS_HPP

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#include "implicit_rkintegrator.hpp"

/**
 * \brief Implements the linear implicit midpoint rule for an autonomous ODE.
 *
 * \tparam Function type for function implementing the rhs function.
 * Must have VectorXd operator()(VectorXd x)
 * \tparam Jacobian type for function implementing the Jacobian of f.
 * Must have MatrixXd operator()(VectorXd x)
 * \param f function handle for rhs in y' = f(y), e.g. implemented
 * using lambda function.
 * \param Jf function handle for Jf, e.g. implemented using lambda function
 * \param T final time T
 * \param y0 initial data y(0) = y0 for y' = f(y)
 * \param N number of steps to perform.
 * \return std::vector<Eigen::VectorXd> containing all steps
 */
/* SAM_LISTING_BEGIN_2 */
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> solve_lin_mid(Function &&f, Jacobian &&Jf,
                                           double T, const Eigen::VectorXd &y0,
                                           unsigned int N) {
  std::vector<Eigen::VectorXd> res;
  // TODO: (12-2.h) Implement the linear implicit mid-point method for
  // an autonomous ODE y' = f(y), y(0) = y0. Return the vector containing
  // all steps including initial and final value.
  // START

  // END
  return res;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void tab_crossprod() {
  // set data
  constexpr double T = 10.;
  constexpr unsigned int N = 128;
  constexpr double c = 1.;
  Eigen::Vector3d y0;
  y0 << 1., 1., 1.;
  Eigen::Vector3d a;
  a << 1., 0., 0.;
  // define rhs
  auto f = [&a, c](const Eigen::Vector3d &y) -> Eigen::Vector3d {
    return a.cross(y) + c * y.cross(a.cross(y));
  };
  // TODO: (12-2.f) solve the cross-product ODE with the implicit RK method
  // defined in solve_imp_mid. Tabulate the norms of the results at all steps.
  // START

  // END
  /* SAM_LISTING_END_3 */

  /* SAM_LISTING_BEGIN_4 */
  // TODO: (12-2.h) solve the cross-product ODE with the implicit RK method
  // defined in solve_lin_mid. Tabulate the norms of the results at all steps.
  // START

  // END
  /* SAM_LISTING_END_4 */
}

#endif
