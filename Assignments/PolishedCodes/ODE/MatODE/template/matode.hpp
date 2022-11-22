#ifndef MATODE_HPP
#define MATODE_HPP

#include <Eigen/Dense>
#include <Eigen/QR>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * @brief One step of the explicit Euler method
 *
 * @param A matrix A of the ODE
 * @param Y0 Initial state
 * @param h step size
 * @return Eigen::MatrixXd the next state
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd eeulstep(const Eigen::MatrixXd &A, const Eigen::MatrixXd &Y0,
                         double h) {
  Eigen::MatrixXd Y1 = Y0;
  // TODO: (11-1.b) Implement ONE step of explicit euler applied to Y0,
  // for the ODE Y' = A*Y
  // START

  // END
  return Y1;
}
/* SAM_LISTING_END_3 */

/**
 * @brief One step of the implicit Euler method
 *
 * @param A matrix A of the ODE
 * @param Y0 Initial state
 * @param h step size
 * @return Eigen::MatrixXd the next state
 */
/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXd ieulstep(const Eigen::MatrixXd &A, const Eigen::MatrixXd &Y0,
                         double h) {
  const unsigned int n = A.rows();
  Eigen::MatrixXd Y1 = Y0;
  // TODO: (11-1.b) Implement ONE step of implicit euler applied to Y0,
  // for the ODE Y' = A*Y
  // START

  // END
  return Y1;
}
/* SAM_LISTING_END_4 */

/**
 * @brief One step of the implicit mid-point method
 *
 * @param A matrix A of the ODE
 * @param Y0 Initial state
 * @param h step size
 * @return Eigen::MatrixXd the next state
 */
/* SAM_LISTING_BEGIN_5 */
Eigen::MatrixXd impstep(const Eigen::MatrixXd &A, const Eigen::MatrixXd &Y0,
                        double h) {
  const unsigned int n = A.rows();
  Eigen::MatrixXd Y1 = Y0;
  // TODO: (11-1.b) Implement ONE step of implicit midpoint rule applied to Y0,
  // for the ODE Y' = A*Y
  // START

  // END
  return Y1;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
std::tuple<double, double, double> checkOrthogonality() {
  constexpr unsigned int n = 3;
  constexpr double h = 0.01;
  Eigen::MatrixXd M(n, n);
  M << 8, 1, 6, 3, 5, 7, 9, 9, 2;

  // Build A
  Eigen::MatrixXd A(n, n);
  A << 0, 1, 1, -1, 0, 1, -1, -1, 0;

  std::vector<double> norms(3, -1.);  //< the output

  // TODO: (11-1.c) Compute and tabulate the Frobenius norms of Y_k'*Y_k - I
  // for 20 steps of eeulstep, ieulstep and impstep.
  // Return the values corresponding to the last step.
  // START

  // END
  return std::make_tuple(norms[0], norms[1], norms[2]);
}
/* SAM_LISTING_END_6 */

#endif
