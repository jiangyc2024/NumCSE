#ifndef MATODE_HPP
#define MATODE_HPP

#include <Eigen/Dense>
#include <Eigen/QR>
#include <iomanip>
#include <iostream>
#include <vector>

//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd eeulstep(const Eigen::MatrixXd &A, const Eigen::MatrixXd &Y0,
                         double h) {
  // TO DO: Implement ONE step of explicit euler applied to Y0,
  // for the ODE Y' = A*Y
  // START
  return Y0;
  // END
}
/* SAM_LISTING_END_3 */

//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXd ieulstep(const Eigen::MatrixXd &A, const Eigen::MatrixXd &Y0,
                         double h) {
  unsigned n = A.rows();
  // TO DO: Implement ONE step of implicit euler applied to Y0,
  // for the ODE Y' = A*Y
  // START
  return Y0;
  // END
}
/* SAM_LISTING_END_4 */

//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
/* SAM_LISTING_BEGIN_5 */
Eigen::MatrixXd impstep(const Eigen::MatrixXd &A, const Eigen::MatrixXd &Y0,
                        double h) {
  unsigned n = A.rows();
  // TO DO: Implement ONE step of implicit midpoint rule applied to Y0,
  // for the ODE Y' = A*Y
  // START
  return Y0;
  // END
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
std::tuple<double, double, double> checkOrthogonality(void) {
  // TO DO (11-1.c): compute and tabulate the Frobenius norms of Y_k'*Y_k - I
  // for 20 steps of eeulstep, ieulstep and impstep.
  // Return the values corresponding to the last step.
  // START
  return std::make_tuple(1.0, 1.0, 1.0);
  // END
}
/* SAM_LISTING_END_6 */

#endif
