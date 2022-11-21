#ifndef ORDNOTALL_H_
#define ORDNOTALL_H_

#include <Eigen/Core>
#include <iomanip>
#include <iostream>
#include <vector>

#include "polyfit.hpp"
#include "rkintegrator.hpp"

/*!
 * \brief errors Approximates the order of convergence of a scheme.
 * The scheme is a Runge Kutta scheme defined by A and b when applied
 * to the first order system y'=f(y), y(0)=y0.
 * The function ouputs error of the solutions at the time T.
 * \tparam Function Type for r.h.s function f.
 * \param f The r.h.s function for the ODE.
 * \param T Final time.
 * \param y0 Initial data.
 * \param A Butcher matrix $A$.
 * \param b Butcher vector $b$.
 */
/* SAM_LISTING_BEGIN_1 */
template <class Function>
double testCvgRKSSM(const Function &f, double T, double y0,
                    const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
  // Helper object carrying out the actual explicit RK-SSM
  RKIntegrator<double> rk(A, b);
  double conv_rate = 0;
  // TO DO: 11-6.a
  // START

  // Vector for collecting errors
  Eigen::VectorXd error(12);
  // Vector for storing No. of timesteps
  Eigen::VectorXd timesteps(12);

  double sum = 0;
  int count = 0;
  bool test = true;

  // Reference numerical solution obtained with $2^{15}$ timesteps
  std::vector<double> y_exact = rk.solve(f, T, y0, std::pow(2, 15));

  for (int k = 0; k < 12; k++) {
    // Number of timesteps
    int M = std::pow(2, k + 1);
    timesteps(k) = M;
    // Solve IVP
    std::vector<double> y1 = rk.solve(f, T, y0, M);
    // Error at final time
    error(k) = std::abs(y1[M] - y_exact[std::pow(2, 15)]);

    std::cout << std::left << std::setfill(' ') << std::setw(3)
              << "M = " << std::setw(7) << M << std::setw(8)
              << "Error = " << std::setw(13) << error[k];

    if (error[k] < 1e-14) {
      test = false;
    }
    if (k > 0 && test) {
      double rk = std::log2(error(k - 1) / error(k));
      std::cout << std::left << std::setfill(' ') << std::setw(10)
                << "Approximated order = " << rk << std::endl;
      sum += rk;
      ++count;
    } else
      std::cout << std::endl;
  }
  std::cout << "Average approximated order = " << sum / count << std::endl
            << std::endl;
  // Compute convergence rate by linear regression
  Eigen::VectorXd fit_coeffs =
      polyfit(timesteps.segment(0, count + 1).array().log(),
              error.segment(0, count + 1).array().log(), 1);
  conv_rate = -fit_coeffs(0);
  // END
  return conv_rate;
}
/* SAM_LISTING_END_1 */

/*!
 * \brief This function compares the convergence rates of four RK single step
 * methods: explicit Euler, trapezoidal rule, RK order 3 and classical RK
 * order 4. Comparison is done for two ODEs: 1. ODE y' = (1-y)y, y(0)=.5 and 2.
 * ODE y' = |1.1 - y| + 1, y(0)=1.
 */
/* SAM_LISTING_BEGIN_2 */
void cmpCvgRKSSM() {
  // TO DO: 11-6.c
  // START
  // Construct data for Butcher schemes
  Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(1, 1);
  Eigen::VectorXd b1(1);
  b1 << 1;

  Eigen::MatrixXd A2 = Eigen::MatrixXd::Zero(2, 2);
  A2(1, 0) = 1;
  Eigen::VectorXd b2(2);
  b2 << .5, .5;

  Eigen::MatrixXd A3 = Eigen::MatrixXd::Zero(3, 3);
  A3(1, 0) = .5;
  A3(2, 0) = -1;
  A3(2, 1) = 2;
  Eigen::VectorXd b3(3);
  b3 << 1. / 6, 2. / 3, 1. / 6;

  Eigen::MatrixXd A4 = Eigen::MatrixXd::Zero(4, 4);
  A4(1, 0) = .5;
  A4(2, 1) = .5;
  A4(3, 2) = 1;
  Eigen::VectorXd b4(4);
  b4 << 1. / 6, 1. / 3, 1. / 3, 1. / 6;

  // First ODE
  std::cout << std::endl
            << "1. ODE y' = (1-y)y, y(0)=.5" << std::endl
            << std::endl;
  double T = 1;
  auto f = [](double y) { return (1. - y) * y; };
  double y0 = .5;
  double rate;

  std::cout << "Explicit Euler" << std::endl;
  std::cout << std::string(50, '_') << std::endl;
  rate = testCvgRKSSM(f, T, y0, A1, b1);
  std::cout << "Estimated rate using linear regression = " << rate << std::endl
            << std::endl;
  std::cout << "Trapezoidal rule" << std::endl;
  std::cout << std::string(50, '_') << std::endl;
  rate = testCvgRKSSM(f, T, y0, A2, b2);
  std::cout << "Estimated rate using linear regression = " << rate << std::endl
            << std::endl;
  std::cout << "RK order 3" << std::endl << std::endl;
  std::cout << std::string(50, '_') << std::endl;
  rate = testCvgRKSSM(f, T, y0, A3, b3);
  std::cout << "Estimated rate using linear regression = " << rate << std::endl
            << std::endl;
  std::cout << "Classical RK order 4" << std::endl << std::endl;
  std::cout << std::string(50, '_') << std::endl;
  rate = testCvgRKSSM(f, T, y0, A4, b4);
  std::cout << "Estimated rate using linear regression = " << rate << std::endl
            << std::endl;

  // Second ODE
  std::cout << std::endl
            << "2. ODE y' = |1.1 - y| + 1, y(0)=1" << std::endl
            << std::endl;
  auto f2 = [](double y) { return std::abs(1.1 - y) + 1.; };
  y0 = 1;

  std::cout << "Explicit Euler" << std::endl << std::endl;
  std::cout << std::string(50, '_') << std::endl;
  rate = testCvgRKSSM(f2, T, y0, A1, b1);
  std::cout << "Estimated rate using linear regression = " << rate << std::endl
            << std::endl;
  std::cout << "Trapezoidal rule" << std::endl << std::endl;
  std::cout << std::string(50, '_') << std::endl;
  rate = testCvgRKSSM(f2, T, y0, A2, b2);
  std::cout << "Estimated rate using linear regression = " << rate << std::endl
            << std::endl;
  std::cout << "RK order 3" << std::endl << std::endl;
  std::cout << std::string(50, '_') << std::endl;
  rate = testCvgRKSSM(f2, T, y0, A3, b3);
  std::cout << "Estimated rate using linear regression = " << rate << std::endl
            << std::endl;
  std::cout << "Classical RK order 4" << std::endl << std::endl;
  std::cout << std::string(50, '_') << std::endl;
  rate = testCvgRKSSM(f2, T, y0, A4, b4);
  std::cout << "Estimated rate using linear regression = " << rate << std::endl
            << std::endl;
  // END
}
/* SAM_LISTING_END_2 */

#endif
