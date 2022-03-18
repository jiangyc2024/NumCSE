#ifndef QUASILIN_HPP
#define QUASILIN_HPP

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <iostream>

/**
 * @brief Computes one step of the fixed point method.
 *
 * @param xk previous step $x^{(k)}$
 * @param b rhs vector $b \in \mathbf{R}^n$
 * @return Eigen::VectorXd next step $x^{(k+1)}$
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd fixed_point_step(const Eigen::VectorXd &xk,
                                 const Eigen::VectorXd &b) {
  Eigen::VectorXd x_new = Eigen::VectorXd::Zero(b.size());
  // TODO: (8-10.b) implement one step of the fixed point iteration
  // START

  // END
  return x_new;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solves the quasilinear system
 *
 * @param rtol relative error tolerance
 * @param atol absolute error tolerance
 * @param b rhs vector $b \in \mathbf{R}^n$
 * @return Eigen::VectorXd next step $x^{(k+1)}$
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solveQuasiLinSystem(double rtol, double atol,
                                    Eigen::VectorXd &b) {
  const int n = b.size();

  Eigen::VectorXd x = b;
  // TODO: (8-10.c) apply the iterative method fixed_point_step to solve
  // the quasi-linear system
  // START

  // END
  return x;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Computes one step of the Newton method
 *
 * @param x previous step $x^{(k)}$
 * @param b rhs vector $b \in \mathbf{R}^n$
 * @return Eigen::VectorXd next step in Newton iteration $x^{(k+1)}$
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd newton_step(const Eigen::VectorXd &x,
                            const Eigen::VectorXd &b) {
  Eigen::VectorXd x_new = Eigen::VectorXd::Zero(b.size());
  // TODO: (8-10.f) implement one Newton step exploiting the
  // Sherman-Morrison-Woodbury formula.
  // START

  // END
  return x_new;
}
/* SAM_LISTING_END_3 */

/**
 * @brief Solves the quasilinear system using Newton's method
 *
 * @param rtol relative error tolerance
 * @param atol absolute error tolerance
 * @param b rhs vector $b \in \mathbf{R}^n$
 * @return Eigen::VectorXd next step $x^{(k+1)}$
 */
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd solveQLSystem_Newton(double rtol, double atol,
                                     Eigen::VectorXd &b) {
  const unsigned int n = b.size();

  Eigen::VectorXd x = b;
  // TODO: (8-10.c) apply the iterative method newton_step to solve
  // the quasi-linear system
  // START

  // END
  return x;
}
/* SAM_LISTING_END_4 */

#endif
