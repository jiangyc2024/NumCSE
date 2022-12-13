#ifndef ODESOLVE_HPP
#define ODESOLVE_HPP

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "matplotlibcpp.h"
#include "polyfit.hpp"

namespace plt = matplotlibcpp;

/**
 * \brief Applies the discrete evolution operator Psi tilde to y0.
 *
 * \tparam DiscEvlOp type for evolution operator (e.g. lambda function type)
 * \param Psi original evolution operator, must have operator(double, const
 * Eigen::VectorXd&)
 * \param p parameter p for construction of Psi tilde
 * \param h step-size
 * \param y0 previous step
 * \return Eigen::VectorXd next state
 */
/* SAM_LISTING_BEGIN_0 */
template <class DiscEvlOp>
Eigen::VectorXd psitilde(DiscEvlOp &&Psi, unsigned int p, double h,
                         const Eigen::VectorXd &y0) {
  Eigen::VectorXd y1 = y0;  // overwrite this
  // TODO: (11-3.b) apply the evolution operator \tilde{\Psi}
  // with step-size h to the value y0
  // START

  // END
  return y1;
}
/* SAM_LISTING_END_0 */

/**
 * \brief Applies psi until final time T.
 *
 * \tparam DiscEvlOp type for evolution operator (e.g. lambda function type)
 * \param Psi original evolution operator, must have operator(double, const
 * Eigen::VectorXd&)
 * \param T final time
 * \param y0 initial data
 * \param N number of steps
 * \return std::vector<Eigen::VectorXd> of states
 */
/* SAM_LISTING_BEGIN_1 */
template <class DiscEvlOp>
std::vector<Eigen::VectorXd> odeintequi(DiscEvlOp &&Psi, double T,
                                        const Eigen::VectorXd &y0,
                                        unsigned int N) {
  std::vector<Eigen::VectorXd> Y;
  // TODO: (11-3.c) Compute y from time 0 to T using N equidistant time steps
  // return a std::vector containing all steps y_0,...,y_N
  // START

  // END
  return Y;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double testcvpExtrapolatedEuler() {
  double conv_rate = 0.;
  constexpr double T = 1.;
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  auto f = [](const Eigen::VectorXd &y) -> Eigen::VectorXd {
    return Eigen::VectorXd::Ones(1) + y * y;
  };
  // TODO: (11-3.d) tabulate the values of the error corresponding to
  // \tilde{\psi}, where \psi is the explicit Euler method.
  // return the empirical convergence rate using polyfit.
  // Hint: first define a lambda for \psi. Then use psitilde to obtain a
  // suitable input for odeintequi.

  // START

  // END
  return conv_rate;
}
/* SAM_LISTING_END_2 */

/**
 * \brief Adaptive timestepper based on psi tilde.
 *
 * \tparam DiscEvlOp type for evolution operator (e.g. lambda function type)
 * \param Psi low-order evolution operator, must have operator(double,
 * const Eigen::VectorXd&)
 * \param T final time
 * \param y0 initial data
 * \param h0 initial step size
 * \param p parameter p for construction of Psi tilde
 * \param reltol relative tolerance for error control
 * \param abstol absolute tolerance for error control
 * \param hmin minimal step size
 * \return std::pair<std::vector<double>, std::vector<Eigen::VectorXd>> of time
 * and corresponding state
 */
/* SAM_LISTING_BEGIN_3 */
template <class DiscEvlOp>
std::pair<std::vector<double>, std::vector<Eigen::VectorXd>> odeintssctrl(
    DiscEvlOp &&Psi, double T, const Eigen::VectorXd &y0, double h0,
    unsigned int p, double reltol, double abstol, double hmin) {
  std::vector<double> t;
  std::vector<Eigen::VectorXd> Y;
  // TODO: (11-3.e)  Compute y from time 0 to T with adaptive time step.
  // Display a warning if the tolerance cannot be met with minimum
  // step size hmin. return a pair of vectors containing the times and
  // the computed values.
  // START

  // END
  return std::make_pair(t, Y);
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
void solveTangentIVP() {
  auto f = [](const Eigen::VectorXd &y) -> Eigen::VectorXd {
    return Eigen::VectorXd::Ones(1) + y * y;
  };
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  plt::figure();
  // TODO: (11-3.f) run the adaptive integration algorithm and plot the
  // resulting values of y(t).
  // Hint: you might use a loop to convert a std::vector<Eigen::VectorXd> into a
  // std::vector<double>, since each Eigen::VectorXd has size 1

  // START

  // END
  plt::savefig("./cx_out/tangent.png");
}
/* SAM_LISTING_END_4 */
#endif
