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

/**
 * @brief Creates the matrix A.
 *
 * @param t $n$ size vector
 * @return Eigen::MatrixXd $n \times 4$ matrix
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd make_A(const Eigen::VectorXd &t) {
  std::size_t n = t.size();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, 4);
  // TODO: (3-3.a) Build the matrix A
  // Hint: Evaluate the functions \phi_j at the time points defined in b

  // START
  for (std::size_t i = 0; i < n; i++) {
    const double t_i = t[i];
    // Initialize row i of the matrix
    A(i, 0) = 1.0 / t_i;
    A(i, 1) = 1.0 / (t_i * t_i);
    A(i, 2) = std::exp(-(t_i - 1.0));
    A(i, 3) = std::exp(-2.0 * (t_i - 1.0));
  }
  // END

  return A;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Computes the least squares estimate for the coefficients $\gamma_i$
 * using the normal equations.
 *
 * @param t $n$ size vector
 * @param f right hand side vector
 * @return Eigen::VectorXd vector $\gamma$ of size 4
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd data_fit_normal(const Eigen::VectorXd &t,
                                const Eigen::VectorXd &f) {
  Eigen::VectorXd gamma = Eigen::VectorXd::Zero(4);

  // TODO: (3-3.a) Solve normal equations to find the coefficients of the
  // linear fitting
  // START
  // Create matrix A
  Eigen::MatrixXd A = make_A(t);

  // Transpose matrix A
  Eigen::MatrixXd At = A.transpose();

  // Compute LHS (A^T*A) of the normal equation
  Eigen::MatrixXd AtA = At * A;

  // Compute RHS (A^T*b) of the normal equation
  Eigen::VectorXd Atf = At * f;

  // Solve normal equation to get coefficients
  // we use the Cholesky decomposition as the matrix AtA is positive
  // semidefinite
  gamma = AtA.ldlt().solve(Atf);
  // END

  return gamma;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Computes the least squares estimate for the coefficients $\gamma_i$
 * using orthogonal transformation techniques.
 *
 * @param t $n$ size vector
 * @param f right hand side vector
 * @return Eigen::VectorXd vector $\gamma$ of size 4
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd data_fit_qr(const Eigen::VectorXd &t,
                            const Eigen::VectorXd &f) {
  Eigen::VectorXd gamma = Eigen::VectorXd::Zero(4);

  // TODO: (3-3.b) Find the coefficients for the linear
  // fitting by means of the QR decomposition of A
  // START
  // Make matrix A
  Eigen::MatrixXd A = make_A(t);

  // Compute QR decomposition
  gamma = A.colPivHouseholderQr().solve(f);
  // END

  return gamma;
}
/* SAM_LISTING_END_3 */

/**
 * @brief Computes f using the estimated parameters $\gamma$ at t_vec.
 *
 * @param gamma size $4$ vector of coefficients
 * @param t_vec vector of evaluation points
 * @return Eigen::VectorXd evaluated function
 */
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd fitted_function(const Eigen::VectorXd &gamma,
                                const Eigen::VectorXd &t_vec) {
  const std::size_t n = t_vec.size();
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n);

  // TODO: (3-3.c) Compute the values of f(t) using the estimated parameters.
  // START
  Eigen::MatrixXd A = make_A(t_vec);
  y = A * gamma;
  // END

  return y;
}
/* SAM_LISTING_END_4 */

/**
 * @brief Computes the squared errors at ten sample points.
 *
 * @param gamma size $4$ vector of coefficients
 * @return Eigen::VectorXd the 10 squared errors
 */
/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd fitting_error(const Eigen::VectorXd &gamma) {
  Eigen::VectorXd err = 1000. * Eigen::VectorXd::Ones(10);
  Eigen::VectorXd t_vec = Eigen::VectorXd::LinSpaced(10, 0.1, 1.0);
  Eigen::VectorXd f(10);
  f << 100., 34., 17., 12., 9., 6., 5., 4., 4., 2.;

  // TODO: (3-3.d) Compute the vector of squared errors of your
  // fit at the provided sample points
  // START
  err = fitted_function(gamma, t_vec) - f;
  err = err.cwiseProduct(err);
  // END

  return err;
}
/* SAM_LISTING_END_5 */

#endif
