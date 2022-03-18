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
  // Define matrix by fixed point $x^k$
  const unsigned int n = xk.size();
  const double nrm = xk.norm();
  Eigen::SparseMatrix<double> A(n, n);
  A.reserve(3 * n);
  for (unsigned int i = 0; i < n; ++i) {
    if (i > 0) A.insert(i, i - 1) = 1;
    A.insert(i, i) = 3 + nrm;
    if (i < n - 1) A.insert(i, i + 1) = 1;
  }
  // solve linear system
  Eigen::SparseLU<Eigen::SparseMatrix<double>> Ax_lu(A);
  x_new = Ax_lu.solve(b);
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
  Eigen::VectorXd x_new(n);
  constexpr unsigned int max_itr = 100;
  for (unsigned int itr = 0; itr < max_itr; ++itr) {
    x_new = fixed_point_step(x, b);

    const double err = (x - x_new).norm();
    x = x_new;

    // Check if Fixed point method is successful
    if (err < atol || err < rtol * x_new.norm()) {
      break;
    }
  }
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
  // Define matrix by fixed point $x^k$
  const unsigned int n = x.size();
  const double nrm = x.norm();
  Eigen::SparseMatrix<double> A(n, n);
  A.reserve(3 * n);
  for (int i = 0; i < n; ++i) {
    if (i > 0) A.insert(i, i - 1) = 1;
    A.insert(i, i) = 3 + nrm;
    if (i < n - 1) A.insert(i, i + 1) = 1;
  }

  // Reuse LU decomposition with SMW
  Eigen::SparseLU<Eigen::SparseMatrix<double>> Ax_lu(A);

  // Solve a bunch of systems
  Eigen::VectorXd Axinv_b = Ax_lu.solve(b);
  Eigen::VectorXd Axinv_x = Ax_lu.solve(x);
  // Next step
  x_new = Axinv_b + Axinv_x * (x.transpose() * (x - Axinv_b)) /
                        (x.norm() + x.dot(Axinv_x));
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
  Eigen::VectorXd x_new(n);
  constexpr unsigned int max_itr = 100;
  for (unsigned int itr = 0; itr < max_itr; ++itr) {
    x_new = newton_step(x, b);

    const double err = (x - x_new).norm();
    x = x_new;

    // Check if Newton method is successful
    if (err < atol || err < rtol * x_new.norm()) {
      break;
    }
  }
  // END
  return x;
}
/* SAM_LISTING_END_4 */

#endif
