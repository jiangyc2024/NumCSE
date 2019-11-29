#ifndef QUASILIN_HPP
#define QUASILIN_HPP

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <iostream>

//! \param[in] b rhs vector $b \in \mathbf{R}^n$
//! \param[in] xk previous step $x^{(k)}$
//! \param[out] x_new next step $x^{(k+1)}$
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd fixed_point_step(const Eigen::VectorXd &xk,
                                 const Eigen::VectorXd &b) {
  Eigen::VectorXd x_new;
  // TO DO (9-10.b): implement one step of the fixed point iteration
  // START
  // Define matrix by fixed point $x^k$
  const int n = xk.size();
  double nrm = xk.norm();
  Eigen::SparseMatrix<double> A(n, n);
  A.reserve(3 * n);
  for (int i = 0; i < n; ++i) {
    if (i > 0)
      A.insert(i, i - 1) = 1;
    A.insert(i, i) = 3 + nrm;
    if (i < n - 1)
      A.insert(i, i + 1) = 1;
  }
  // solve linear system
  Eigen::SparseLU<Eigen::SparseMatrix<double>> Ax_lu(A);
  x_new = Ax_lu.solve(b);
  // END
  return x_new;
}
/* SAM_LISTING_END_1 */

//! \param[in] rtol relative error tolerance
//! \param[in] atol absolute error tolerance
//! \param[in] b rhs vector $b \in \mathbf{R}^n$
//! \param[out] x next step $x^{(k+1)}$
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solveQuasiLinSystem(double rtol, double atol,
                                    Eigen::VectorXd &b) {
  int n = b.size();

  Eigen::VectorXd x(n);
  // TO DO (9-10.c): apply the iterative method fixed_point_step to solve
  // the quasi-linear system
  // START
  x = b;
  Eigen::VectorXd x_new(n);
  int max_itr = 100;
  for (int itr = 0; itr < max_itr; ++itr) {

    x_new = fixed_point_step(x, b);

    double err = (x - x_new).norm();
    x = x_new;

    if (err < atol || err < rtol * x_new.norm()) {
      std::cout << "Fixed point method successful" << std::endl;
      break;
    }
  }
  // END
  return x;
}
/* SAM_LISTING_END_2 */

//! \param[in] b rhs vector $b \in \mathbf{R}^n$
//! \param[in] x previous step $x^{(k)}$
//! \param[out] x_new next step in Newton iteration $x^{(k+1)}$
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd newton_step(const Eigen::VectorXd &x,
                            const Eigen::VectorXd &b) {
  Eigen::VectorXd x_new;
  // TO DO (9-10.f): implement one Newton step exploiting the
  // Sherman-Morrison-Woodbury formula.
  // START
  // Define matrix by fixed point $x^k$
  const int n = x.size();
  double nrm = x.norm();
  Eigen::SparseMatrix<double> A(n, n);
  A.reserve(3 * n);
  for (int i = 0; i < n; ++i) {
    if (i > 0)
      A.insert(i, i - 1) = 1;
    A.insert(i, i) = 3 + nrm;
    if (i < n - 1)
      A.insert(i, i + 1) = 1;
  }

  // Reuse LU decomposition with SMW
  Eigen::SparseLU<Eigen::SparseMatrix<double>> Ax_lu(A);

  // Solve a bunch of systems
  Eigen::VectorXd Axinv_b = Ax_lu.solve(b);
  Eigen::VectorXd Axinv_x = Ax_lu.solve(x);
  // Next step
  x_new = Axinv_b + Axinv_x * ( x.transpose() * (x - Axinv_b)) /
                        (x.norm() + x.dot(Axinv_x));
  // END
  return x_new;
}
/* SAM_LISTING_END_3 */

//! \param[in] rtol relative error tolerance
//! \param[in] atol absolute error tolerance
//! \param[in] b rhs vector $b \in \mathbf{R}^n$
//! \param[out] x next step $x^{(k+1)}$
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd solveQLSystem_Newton(double rtol, double atol,
                                     Eigen::VectorXd &b) {
  int n = b.size();

  Eigen::VectorXd x(n);
  // TO DO (9-10.c): apply the iterative method newton_step to solve
  // the quasi-linear system
  // START
  x = b;
  Eigen::VectorXd x_new(n);
  int max_itr = 100;
  for (int itr = 0; itr < max_itr; ++itr) {

    x_new = newton_step(x, b);

    double err = (x - x_new).norm();
    x = x_new;

    if (err < atol || err < rtol * x_new.norm()) {
      std::cout << "Newton method successful" << std::endl;
      break;
    }
  }
  // END
  return x;
}
/* SAM_LISTING_END_4 */

#endif
