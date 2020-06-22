#ifndef QUASILIN_HPP
#define QUASILIN_HPP

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <iostream>

//! \param[in] b rhs vector $b \in \mathbf{R}^n$
//! \param[in] xk previous step $x^{(k)}$
//! \param[out] x_new next step $x^{(k+1)}$
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd fixed_point_step(const Eigen::VectorXd &xk, const Eigen::VectorXd &b) {
  Eigen::VectorXd x_new;
  // TO DO (9-10.b): implement one step of the fixed point iteration
  // START
  
  // END
  return x_new;
}
/* SAM_LISTING_END_1 */

//! \param[in] rtol relative error tolerance
//! \param[in] atol absolute error tolerance
//! \param[in] b rhs vector $b \in \mathbf{R}^n$
//! \param[out] x next step $x^{(k+1)}$
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solveQuasiLinSystem(double rtol, double atol, Eigen::VectorXd &b) {
  int n = b.size();

  Eigen::VectorXd x(n);
  // TO DO (9-10.c): apply the iterative method fixed_point_step to solve
  // the quasi-linear system
  // START
  
  // END
  return x;
}
/* SAM_LISTING_END_2 */

//! \param[in] b rhs vector $b \in \mathbf{R}^n$
//! \param[in] x previous step $x^{(k)}$
//! \param[out] x_new next step in Newton iteration $x^{(k+1)}$
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd newton_step(const Eigen::VectorXd &x, const Eigen::VectorXd &b) {
  Eigen::VectorXd x_new;
  // TO DO (9-10.f): implement one Newton step exploiting the
  // Sherman-Morrison-Woodbury formula.
  // START
  
  // END
  return x_new;
}
/* SAM_LISTING_END_3 */

//! \param[in] rtol relative error tolerance
//! \param[in] atol absolute error tolerance
//! \param[in] b rhs vector $b \in \mathbf{R}^n$
//! \param[out] x next step $x^{(k+1)}$
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd solveQLSystem_Newton(double rtol, double atol, Eigen::VectorXd &b) {
  int n = b.size();

  Eigen::VectorXd x(n);
  // TO DO (9-10.c): apply the iterative method newton_step to solve
  // the quasi-linear system
  // START
  
  // END
  return x;
}
/* SAM_LISTING_END_4 */

#endif
