#ifndef GRADIENTFLOW_HPP
#define GRADIENTFLOW_HPP

#include <Eigen/Dense>
#include <array>
#include <iostream>
#include <vector>

/**
 * \brief Implements the Butcher matrix of the SDIRK method
 *
 * \return Eigen::MatrixXd Butcher matrix
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd ButcherMatrix() {
  Eigen::MatrixXd A(6, 5);
  // clang-format off
  A <<       0.25,          0.,       0.,       0.,   0.,
    0.5,        0.25,       0.,       0.,   0.,
    17./50.,     -1./25.,     0.25,       0.,   0.,
    371./1360., -137./2720., 15./544.,     0.25,   0.,
    25./24.,    -49./48., 125./16., -85./12., 0.25,
    25./24.,    -49./48., 125./16., -85./12., 0.25;
  // clang-format on
  return A;
}
/* SAM_LISTING_END_0 */

/**
 * \brief Solves one stage equation using the (undamped) Newton method
 *
 * \tparam Functor function handle
 * \tparam Jacobian function handle
 * \param f r.h.s. function
 * \param J Jacobian of f
 * \param y initial state
 * \param b r.h.s.
 * \param h step size
 * \param rtol relative tolerance
 * \param atol absolute tolerance
 * \return Eigen::VectorXd stage
 */
/* SAM_LISTING_BEGIN_1 */
template <typename Functor, typename Jacobian>
Eigen::VectorXd solveGenStageEquation(Functor &&f, Jacobian &&J,
                                      const Eigen::VectorXd &y,
                                      const Eigen::VectorXd &b, double h,
                                      double rtol = 1E-6, double atol = 1E-8) {
  Eigen::VectorXd g = Eigen::VectorXd::Zero(
      y.size());  // initial guess g=0 for Newton iteration.

  // TODO: (12-8.e; optional) Solve one stage equation using the Newton
  // iteration.
  // START

  // END

  return g;
}
/* SAM_LISTING_END_1 */

/**
 * \brief Computes the five stages of the SDIRK method
 *
 * \tparam Functor function handle
 * \tparam Jacobian function handle
 * \param f r.h.s. function
 * \param J Jacobian of f
 * \param y initial state
 * \param h step size
 * \param rtol relative tolerance
 * \param atol absolute tolerance
 * \return std::array<Eigen::VectorXd, 5> of 5 stages
 */
/* SAM_LISTING_BEGIN_2 */
template <typename Functor, typename Jacobian>
std::array<Eigen::VectorXd, 5> computeStages(Functor &&f, Jacobian &&J,
                                             const Eigen::VectorXd &y, double h,
                                             double rtol = 1E-6,
                                             double atol = 1E-8) {
  std::array<Eigen::VectorXd, 5> G;  // array of stages
  const unsigned int d = y.size();
  Eigen::MatrixXd Coeffs = ButcherMatrix();
  // TODO: (12-8.e) Solve for the stages of the SDIRK method using the
  // (undamped) Newton method.
  // Hint: You may implement solveGenStageEquation as
  // an auxiliary function.
  // START

  // END
  return G;
}
/* SAM_LISTING_END_2 */

/**
 * \brief Perform a discrete evolution using the SDIRK method.
 *
 * \tparam Functor function handle
 * \tparam Jacobian function handle
 * \param f r.h.s. function of the ODE
 * \param J Jacobian of f
 * \param y initial state
 * \param h step size
 * \param rtol relative tolerance
 * \param atol absolute tolerance
 * \return Eigen::VectorXd evolved state
 */
/* SAM_LISTING_BEGIN_5 */
template <typename Functor, typename Jacobian>
Eigen::VectorXd discEvolSDIRK(Functor &&f, Jacobian &&J,
                              const Eigen::VectorXd &y, double h,
                              double rtol = 1E-6, double atol = 1E-8) {
  // The b weights are in the last row of Coeffs.
  Eigen::MatrixXd Coeffs = ButcherMatrix();
  Eigen::VectorXd Psi;
  // TODO: (12-8.f) Apply the discrete evolution operator to the ODE described
  // by f.
  // START

  // END
  return Psi;
}
/* SAM_LISTING_END_5 */

/**
 * \brief Solves the gradient flow ODE as given in the problem description.
 *
 * \param d vector
 * \param lambda real value
 * \param y initial state
 * \param T final time
 * \param M number of steps to take
 * \return std::vector<Eigen::VectorXd> vector of states
 */
/* SAM_LISTING_BEGIN_6 */
std::vector<Eigen::VectorXd> solveGradientFlow(const Eigen::VectorXd &d,
                                               double lambda,
                                               const Eigen::VectorXd &y,
                                               double T, unsigned int M) {
  std::vector<Eigen::VectorXd> Y(M + 1);
  // TODO: (12-8.i) Solve the Gradient flow ODE with the vector potential given
  // in 12.8.8
  // START

  // END
  return Y;
}
/* SAM_LISTING_END_6 */

#endif