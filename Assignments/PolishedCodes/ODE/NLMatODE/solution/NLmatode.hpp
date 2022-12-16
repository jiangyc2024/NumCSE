#ifndef NLMATODE_HPP
#define NLMATODE_HPP

#include <iomanip>
#include <iostream>

#include "ode45.hpp"
#include "polyfit.hpp"

/**
 * \brief Finds an approximation of the matrix IVP $Y' = -(Y-Y')*Y$ at time $T$
 *
 * \param Y0 Initial data Y(0) (as matrix)
 * \param T final time of simulation
 * \return Eigen::MatrixXd solution at final time
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd matode(const Eigen::MatrixXd &Y0, double T) {
  Eigen::MatrixXd YT = Y0;  // overwrite this
  // TODO: (11-5.a): use the ode45 class to find an approximation
  // of the matrix IVP $Y' = -(Y-Y')*Y$ at time $T$
  // START
  // Define the RHS
  auto F = [](const Eigen::MatrixXd &M) { return -(M - M.transpose()) * M; };
  ode45<Eigen::MatrixXd> O(F);

  // Set tolerances
  O.options.atol = 1e-10;
  O.options.rtol = 1e-8;

  // Return only matrix at $T$, (solution is vector
  // of pairs $(y(t_k), t_k)$ for each step k
  YT = O.solve(Y0, T).back().first;
  // END
  return YT;
}
/* SAM_LISTING_END_1 */

/**
 * \brief Checks if invariant $Y'*Y$ is preserved
 *
 * \param M Initial data Y(0) (as matrix)
 * \param T final time of simulation
 * \return true if invariant $Y' * Y$ is preserved
 * \return false otherwise
 */
/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const Eigen::MatrixXd &M, double T) {
  bool inv_holds = false;
  // TODO: (11-5.c) check if $Y'*Y$ is preserved at the time $T$ by matode.
  // START
  Eigen::MatrixXd N = matode(M, T);

  if ((N.transpose() * N - M.transpose() * M).norm() <
      10 * std::numeric_limits<double>::epsilon() * M.norm()) {
    inv_holds = true;
  } else {
    inv_holds = false;
  }
  // END
  return inv_holds;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
double cvgDiscreteGradientMethod() {
  double conv_rate = 0;
  // TODO: (11-5.d) compute the fitted convergence rate of the Discrete
  // gradient method. Also tabulate the values M and the errors.
  // START
  constexpr double T = 1.0;
  // initial value
  Eigen::MatrixXd Y0 = Eigen::MatrixXd::Zero(5, 5);
  Y0(4, 0) = 1;
  for (unsigned int i = 0; i < 4; ++i) {
    Y0(i, i + 1) = 1;
  }
  // reference solution
  Eigen::MatrixXd Y_ex = matode(Y0, T);

  // define the rhs
  auto F = [](const Eigen::MatrixXd &M) { return -(M - M.transpose()) * M; };

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(5, 5);
  Eigen::ArrayXd MM(8);
  Eigen::ArrayXd err(8);

  std::cout << "Error for equidistant steps:" << std::endl;
  std::cout << "M"
            << "\t"
            << "Error" << std::endl;
  for (unsigned int i = 0; i < 8; ++i) {
    const unsigned int M = 10 * std::pow(2, i);
    const double h = T / M;
    MM(i) = M;
    Eigen::MatrixXd Y = Y0;
    for (unsigned int j = 0; j < M; ++j) {
      Eigen::MatrixXd Ystar = Y + 0.5 * h * F(Y);

      Eigen::MatrixXd Yinc = 0.5 * h * (Ystar - Ystar.transpose());
      Y = (I + Yinc).lu().solve((I - Yinc) * Y);
    }
    err(i) = (Y - Y_ex).norm();
    std::cout << M << "\t" << err(i) << std::endl;
  }

  // compute fitted rate
  Eigen::VectorXd coeffs = polyfit(MM.log(), err.log(), 1);
  conv_rate = -coeffs(0);
  // END
  return conv_rate;
}
/* SAM_LISTING_END_3 */

#endif
