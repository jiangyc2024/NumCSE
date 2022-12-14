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
  // Initial step size
  const double h = T / N;
  // Store initial data
  res.push_back(y0);

  // Initialize some memory to store temporary values
  Eigen::VectorXd ytemp1 = y0;
  Eigen::VectorXd ytemp2 = y0;
  // Pointers to swap previous value
  Eigen::VectorXd *yold = &ytemp1;
  Eigen::VectorXd *ynew = &ytemp2;
  Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(3, 3);

  // Loop over all fixed steps
  for (unsigned int k = 0; k < N; ++k) {
    // Compute, save and swap next step
    *ynew = *yold + h * (eye - h / 2. * Jf(*yold)).lu().solve(f(*yold));
    res.push_back(*ynew);
    std::swap(yold, ynew);
  }
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
  // define Jacobian of rhs
  auto Jf = [&a](const Eigen::Vector3d &y) {
    Eigen::Matrix3d temp;
    temp << -c * (a(1) * y(1) + a(2) * y(2)),
        c * (2 * a(0) * y(1) - a(1) * y(0)) - a(2),
        a(1) + c * (2 * a(0) * y(2) - a(2) * y(0)),
        a(2) - c * (a(0) * y(1) - 2 * a(1) * y(0)),
        -c * (a(0) * y(0) + a(2) * y(2)),
        c * (2 * a(1) * y(2) - a(2) * y(1)) - a(0),
        -a(1) - c * (a(0) * y(2) - 2 * a(2) * y(0)),
        a(0) - c * (a(1) * y(2) - 2 * a(2) * y(1)),
        -c * (a(0) * y(0) + a(1) * y(1));
    return temp;
  };

  // Solving the ODE using implicit midpoint method
  std::vector<Eigen::VectorXd> res_imp(N + 1);
  // Initialize implicit RK with Butcher scheme
  constexpr unsigned int s = 1;
  // Initialize coefficients for the implicit midpoint method
  Eigen::MatrixXd A(s, s);
  Eigen::VectorXd b(s);
  A << 1. / 2.;
  b << 1.;
  implicitRKIntegrator RK(A, b);
  res_imp = RK.solve(f, Jf, T, y0, N);

  // Tabulating the norm of generated states
  std::cout << "1. Implicit midpoint method" << std::endl;
  std::cout << std::setw(10) << "t" << std::setw(15) << "norm(y(t))"
            << std::endl;

  for (unsigned int i = 0; i < N + 1; ++i) {
    std::cout << std::setw(10) << T * i / N << std::setw(15)
              << res_imp[i].norm() << std::endl;
  }
  // END
  /* SAM_LISTING_END_3 */

  /* SAM_LISTING_BEGIN_4 */
  // TODO: (12-2.h) solve the cross-product ODE with the implicit RK method
  // defined in solve_lin_mid. Tabulate the norms of the results at all steps.
  // START
  std::vector<Eigen::VectorXd> res_lin = solve_lin_mid(f, Jf, T, y0, N);

  std::cout << "\n2. Linear implicit midpoint method" << std::endl;
  std::cout << std::setw(10) << "t" << std::setw(15) << "norm(y(t))"
            << std::endl;
  for (int i = 0; i < N + 1; ++i) {
    std::cout << std::setw(10) << T * i / N << std::setw(15)
              << res_lin[i].norm() << std::endl;
  }
  // END
  /* SAM_LISTING_END_4 */
}

#endif
