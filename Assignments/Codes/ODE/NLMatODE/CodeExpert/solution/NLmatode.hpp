#ifndef NLMATODE_HPP
#define NLMATODE_HPP

#include "ode45.hpp"
#include "polyfit.hpp"
using namespace Eigen;

//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
/* SAM_LISTING_BEGIN_1 */
MatrixXd matode(const MatrixXd &Y0, double T) {
  // TO DO (12-5.a): use the ode45 class to find an approximation
  // of the matrix IVP $Y' = -(Y-Y')*Y$ at time $T$

  // START
  // Define the RHS
  auto F = [](const MatrixXd &M) { return -(M - M.transpose()) * M; };
  ode45<MatrixXd> O(F);

  // Set tolerances
  O.options.atol = 10e-10;
  O.options.rtol = 10e-8;

  // Return only matrix at $T$, (solution is vector
  // of pairs $(y(t_k), t_k)$ for each step k
  MatrixXd YT = O.solve(Y0, T).back().first;
  // END
  return YT;
}
/* SAM_LISTING_END_1 */

//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const MatrixXd &M, double T) {
  // TO DO (12-5.c): check if $Y'*Y$ is preserved at the time $T$ by matode.
  // START
  MatrixXd N = matode(M, T);

  if ((N.transpose() * N - M.transpose() * M).norm() <
      10 * std::numeric_limits<double>::epsilon() * M.norm()) {
    return true;
  } else {
    return false;
  }
  // END
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
double cvgDiscreteGradientMethod(void) {
  // TO DO (12-5.d): compute the fitted convergence rate of the Discrete
  // gradient method. Also tabulate the values M and the errors.
  double conv_rate;
  // START
  double T = 1.0;
  // initial value
  MatrixXd Y0 = MatrixXd::Zero(5, 5);
  Y0(4, 0) = 1;
  for (unsigned int i = 0; i < 4; ++i) {
    Y0(i, i + 1) = 1;
  }
  // reference solution
  MatrixXd Y_ex = matode(Y0, T);

  // define the rhs
  auto F = [](const MatrixXd &M) { return -(M - M.transpose()) * M; };

  MatrixXd I = MatrixXd::Identity(5, 5);
  // Auxiliary arrays for computing the convergence rate
  ArrayXd MM(8);
  ArrayXd err(8);

  std::cout << "Error for equidistant steps:" << std::endl;
  std::cout << "M\tError" << std::endl;
  // Outer loop for studying convergence
  for (unsigned int i = 0; i < 8; ++i) {
    unsigned int M = 10 * std::pow(2, i);
    double h = T / M;
    MM(i) = M;
    MatrixXd Y = Y0;
    // Timestepping loop
    for (unsigned int j = 0; j < M; ++j) {
      // Auxiliary matrix Y_*
      MatrixXd Ystar = Y + 0.5 * h * F(Y);
      // Second auxiliary matrix: h/2*(Y_* - Y_*^T)
      MatrixXd Yinc = 0.5 * h * (Ystar - Ystar.transpose());
      // Recursion of the discrete gradient method
      Y = (I + Yinc).lu().solve((I - Yinc) * Y);
    }
    // Compute error norm at final time
    err(i) = (Y - Y_ex).norm();
    std::cout << M << "\t" << err(i) << std::endl;
  }
  // rate of algebraic convergence by linear regression
  VectorXd coeffs = polyfit(MM.log(), err.log(), 1);
  conv_rate = -coeffs(0);
  // END
  return conv_rate;
}
/* SAM_LISTING_END_3 */

#endif
