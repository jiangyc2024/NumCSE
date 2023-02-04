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

  // Need to solve the equation lhs(g) = g - h*f(y+g)/4 - b = 0.
  auto lhs = [f, y, b, h](const Eigen::VectorXd &g) {
    // lhs(g) = g - h*f(y+g)/4 - b
    Eigen::VectorXd val = g - 0.25 * h * f(y + g) - b;
    return val;
  };
  auto Jlhs = [J, y, h](const Eigen::VectorXd &g) {
    // lhs(g) = g - h*f(y+g)/4 - b, so the Jacobian is
    // Jlhs(g) = Id - h*Jf(y+g)/4
    int dim = y.size();
    Eigen::MatrixXd Jval =
        Eigen::MatrixXd::Identity(dim, dim) - 0.25 * h * J(y + g);
    return Jval;
  };

  // Perform Newton iterations:
  Eigen::VectorXd delta =
      -Jlhs(g).lu().solve(lhs(g));  // Newton correction term.
  int iter = 0,
      maxiter = 100;  // If correction based termination does not work.
  while (delta.norm() > atol && delta.norm() > rtol * g.norm() &&
         iter < maxiter) {
    g = g + delta;
    delta = -Jlhs(g).lu().solve(lhs(g));
    iter++;
  }
  g = g + delta;  // Perform the final step
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
  for (unsigned int i = 0; i < 5; i++) {
    // Calculate coefficients from previous stages:
    Eigen::VectorXd b = Eigen::VectorXd::Zero(d);
    for (unsigned int j = 0; j < i; j++) {
      b += Coeffs(i, j) * f(y + G.at(j));
    }
    b *= h;
    // Compute stage i and store in G.at(i):
    G.at(i) = solveGenStageEquation(f, J, y, b, h, rtol, atol);
  }
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
  const unsigned int nstages = Coeffs.cols();
  Psi = y;

  // Using computeStages()
  std::array<Eigen::VectorXd, 5> G = computeStages(f, J, y, h, rtol, atol);
  for (int i = 0; i < nstages; i++)
    Psi += h * Coeffs(nstages, i) * f(y + G.at(i));
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
  // Define the right hand side of the ODE y' = f(y), and the
  // Jacobian of f.
  auto f = [d, lambda](const Eigen::VectorXd &yt) {
    Eigen::VectorXd val =
        -2. * std::cos(yt.squaredNorm()) * yt - 2. * lambda * d.dot(yt) * d;
    return val;
  };
  auto J = [d, lambda](const Eigen::VectorXd &yt) {
    const unsigned int dim = yt.size();
    Eigen::MatrixXd term1 =
        4. * std::sin(yt.squaredNorm()) * yt * yt.transpose();
    Eigen::MatrixXd term2 =
        -2. * std::cos(yt.squaredNorm()) * Eigen::MatrixXd::Identity(dim, dim);
    Eigen::MatrixXd term3 = -2. * lambda * d * d.transpose();
    Eigen::MatrixXd Jval = term1 + term2 + term3;
    return Jval;
  };

  // Split the interval [0,T] into M intervals of size h.
  const double h = T / M;
  Eigen::VectorXd yt = y;
  Y.at(0) = y;
  // Evolve up to time T:
  for (unsigned int i = 1; i <= M; i++) {
    yt = discEvolSDIRK(f, J, yt, h);
    Y.at(i) = yt;
  }
  // END
  return Y;
}
/* SAM_LISTING_END_6 */

#endif