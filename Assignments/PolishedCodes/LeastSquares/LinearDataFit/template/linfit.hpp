#ifndef LINFIT_HPP
#define LINFIT_HPP

////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

/* @param[in] b $n$ size vector
 * @return A $n \times 4$ matrix
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd make_A(const Eigen::VectorXd &b) {
  size_t n = b.size();
  Eigen::MatrixXd A(n, 4);
  // TODO (3-3.a) Build the matrix A
  // Hint: Evaluate the functions \phi_j at the time points defined in b

  // START

  // END

  return A;
}
/* SAM_LISTING_END_1 */

/* @param[in] b $n$ size vector
 * @param[in] t $n$ size vector
 * @return gamma $4$ size vector
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd data_fit_normal(const Eigen::VectorXd &t_vec,
                                const Eigen::VectorXd &b) {
  // TODO (3-3.a) Solve normal equations to find the coefficients of the
  // linear fitting
  Eigen::VectorXd gamma(4);

  // START
 
  // END

  return gamma;
}
/* SAM_LISTING_END_2 */

/* @param[in] b $n$ size vector
 * @param[in] t $n $ size vector
 * @return gamma $4$ size vector
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd data_fit_qr(const Eigen::VectorXd &t_vec,
                            const Eigen::VectorXd &b) {
  // TODO (3-3.b) Find the coefficients for the linear
  // fitting by means of the QR decomposition of A
  Eigen::VectorXd gamma(4);

  // START
 
  // END

  return gamma;
}
/* SAM_LISTING_END_3 */

/* @param[in] gamma $4$ size column vector
 * @param[in] t vector
 * @return y vector of size = t.size()
 */
// Note: the code will not run until this function is implemented
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd fitted_function(const Eigen::VectorXd &gamma,
                                const Eigen::VectorXd &t_vec) {
  // TODO (3-3.c): Define the data for the first plot by
  // evaluating the function f at the grid defined by t_vec

  unsigned int n = t_vec.size();
  Eigen::VectorXd y(n);

  // START

  // END

  return y;
}
/* SAM_LISTING_END_4 */

/* @param[in] gamma $4$ size column vector
 * @return err $n$ size vector
 */
/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd fitting_error(const Eigen::VectorXd &gamma) {
  // TODO (3-3.c): Compute the vector of squared errors of your
  // fit at the provided sample points
  Eigen::VectorXd err(10);

  // START
 
  // END

  return err;
}
/* SAM_LISTING_END_5 */

#endif
