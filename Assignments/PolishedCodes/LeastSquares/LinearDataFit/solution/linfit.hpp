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
  for (size_t i = 0; i < n; i++) {
    double t_vec = b(i);
    A(i, 0) = 1.0 / t_vec;
    A(i, 1) = 1.0 / (t_vec * t_vec);
    A(i, 2) = std::exp(-(t_vec - 1.0));
    A(i, 3) = std::exp(-2.0 * (t_vec - 1.0));
  }
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
  // Create matrix A
  Eigen::MatrixXd A = make_A(t_vec);

  // Transpose matrix A
  Eigen::MatrixXd At = A.transpose();

  // Compute LHS (A^T*A) of the normal equation
  Eigen::MatrixXd AtA = At * A;

  // Compute RHS (A^T*b) of the normal equation
  Eigen::VectorXd Atb = At * b;

  // Solve normal equation to get coefficients
  gamma = AtA.ldlt().solve(Atb);
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
  // Make matrix A
  Eigen::MatrixXd A = make_A(t_vec);

  // Compute QR decomposition
  gamma = A.colPivHouseholderQr().solve(b);
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

  Eigen::VectorXd y;

  // START
  Eigen::MatrixXd A = make_A(t_vec);
  y = A * gamma;
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
  Eigen::VectorXd err;

  // START
  Eigen::VectorXd t_vec = Eigen::VectorXd::LinSpaced(10, 0.1, 1.0);
  Eigen::VectorXd f(10);
  f << 100., 34., 17., 12., 9., 6., 5., 4., 4., 2.;

  // Compute error
  err = fitted_function(gamma, t_vec) - f;
  err = err.cwiseProduct(err);
  // END

  return err;
}
/* SAM_LISTING_END_5 */

#endif
