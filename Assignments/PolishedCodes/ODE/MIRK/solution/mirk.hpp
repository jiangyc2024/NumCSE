#ifndef MIRK_HPP
#define MIRK_HPP

/**
 * \file mirk.hpp
 * \brief NPDE homework MIRK code
 * \author Unknown, Oliver Rietmann
 * \date 04.04.2021
 * \copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/LU>
#include <iomanip>
#include <iostream>

/**
 * \brief Perform 2 steps of the Newton method applied to F and its Jacobian DF
 *
 * \tparam Func suitable functor for F
 * \tparam Jac suitable functor for DF
 * \param F function F
 * \param DF Jacobian of function F
 * \param z initial guess
 * \return Eigen::VectorXd z after applying two steps of Newton's method
 */
/* SAM_LISTING_BEGIN_0 */
template <class Func, class Jac>
Eigen::VectorXd Newton2Steps(Func &&F, Jac &&DF, Eigen::VectorXd z) {
  // TODO: (12-6.e) Perform two steps of Newton's method.
  // START
  // First Newton step
  z = z - DF(z).lu().solve(F(z));
  // Second Newton step
  z = z - DF(z).lu().solve(F(z));
  // END
  return z;
}
/* SAM_LISTING_END_0 */

/**
 * \brief Perform a single step of the MIRK scheme applied to the scalar ODE
 * y' = f(y)
 *
 * \tparam Func suitable functor for f
 * \tparam Jac suitable functor for Jacobian of f
 * \param f r.h.s. function
 * \param df Jacobian of f
 * \param y0 initial state
 * \param h step size
 * \return double next state
 */
/* SAM_LISTING_BEGIN_1 */
template <class Func, class Jac>
double MIRKStep(Func &&f, Jac &&df, double y0, double h) {
  // Coefficients of MIRK
  constexpr double v1 = 1.0;
  constexpr double v2 = 344.0 / 2025.0;
  constexpr double d21 = -164.0 / 2025.0;
  constexpr double b1 = 37.0 / 82.0;
  constexpr double b2 = 45.0 / 82.0;

  double y1 = y0;  // next state; overwrite this
  // TODO: (12-6.f) Perform one step of the MIRK scheme applied to the scalar
  // ODE given by f.
  // START

  // F derived from MIRK scheme (vector valued)
  auto F = [&f, &y0, &h](Eigen::Vector3d z) -> Eigen::Vector3d {
    Eigen::Vector3d ret;
    ret << z(0) - (1 - v1) * y0 - v1 * z(2),
        z(1) - (1 - v2) * y0 - v2 * z(2) - h * d21 * f(z(0)),
        z(2) - y0 - h * (b1 * f(z(0)) + b2 * f(z(1)));
    return ret;
  };
  // Jacobian of F (matrix-valued)
  auto DF = [&df, &h, &v1, &v2](Eigen::Vector3d z) -> Eigen::Matrix3d {
    Eigen::Matrix3d M;
    M << 0, 0, v1, h * d21 * df(z(0)), 0, v2, h * b1 * df(z(0)),
        h * b2 * df(z(1)), 0;
    return Eigen::Matrix3d::Identity() - M;
  };
  // Initial Guess for Newton steps
  Eigen::Vector3d z = {0.0, 0.0, y0};

  // Approximate z, s.t. F(z) = 0
  z = Newton2Steps(F, DF, z);
  y1 = z(2);
  // END
  return y1;
}
/* SAM_LISTING_END_1 */

/**
 * \brief Solve an ODE y' = f(y) using MIRK scheme on equidistant steps,
 * return the approximation of y(T)
 *
 * \tparam Func suitable functor for f
 * \tparam Jac suitable functor for the Jacobian of f
 * \param f r.h.s. f
 * \param df Jacobian of f
 * \param y0 initial state
 * \param T final time
 * \param M number of steps
 * \return double final state at T
 */
/* SAM_LISTING_BEGIN_2 */
template <class Func, class Jac>
double MIRKSolve(Func &&f, Jac &&df, double y0, double T, unsigned int M) {
  // Will contain next step
  double y = y0;
  // TODO: (12-6.g) Solve the ODE given by f until time T
  // START

  // Step size
  const double h = T / M;
  // Perform N quidistant steps
  for (unsigned int i = 0; i < M; ++i) {
    y = MIRKStep(f, df, y, h);
  }
  // Return final value at t = T
  // END
  return y;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void cvgMIRK() {
  // TODO: (12-6.h) Tabulate the errors obtained from using the MIRK scheme on
  // the given ODE.
  // START

  // r.h.s
  auto f = [](double y) -> double { return 1 + y * y; };
  // Jacobian of $f$
  auto df = [](double y) -> double { return 2 * y; };
  // Initial data
  constexpr double y0 = 0.;
  // Final time
  constexpr double T = 1.;
  // Exact solution at t = T = 1
  const double yex = std::tan(T);

  std::cout << "Convergence MIRK for IVP y' = 1+y^2 " << std::endl;
  // Table header
  std::cout << "M"
            << "\t"
            << "yend"
            << "\t"
            << "err" << std::endl;
  for (unsigned int M = 4; M <= 512; M *= 2) {
    // Solve up to time T = 1, using M equidistant steps
    const double yend = MIRKSolve(f, df, y0, T, M);
    // Compute error
    const double err = std::abs(yex - yend);
    // Print table
    std::cout << M << "\t" << yend << "\t" << err << std::endl;
  }
  // END
}
/* SAM_LISTING_END_3 */

#endif
