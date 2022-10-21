#ifndef GRAMSCHMIDT_HPP
#define GRAMSCHMIDT_HPP
////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>

/**
 * @brief Performs Gram-Schidt orthonormalization
 * Given a matrix $\mathbf{A}$ of linearly independent columns,
 * returns the result of a Gram-Schmidt orthonormalization.
 * Unstable GS algorithm: output is prone to cancellation issues.
 *
 * @param A Matrix of linearly independent columns
 * @return Eigen::MatrixXd with ONB of $span(a_1, \cdots, a_n)$ as columns
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd gram_schmidt(const Eigen::MatrixXd &A) {
  // We create a matrix Q with the same size and data of A
  Eigen::MatrixXd Q(A);

  // The first vector just gets normalized
  Q.col(0).normalize();

  // TODO: (1-2.b) Implement the gram\_schmidt procedure by iterating over all
  // other columns of A.
  // START
  for (unsigned int j = 1; j < A.cols(); ++j) {
    // See eigen documentation for usage of col and leftCols
    Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));

    // Normalize vector, if possible
    // (otherwise it means columns of $\mathbf{A}$ are
    // almost linearly dependent)
    constexpr double eps = std::numeric_limits<double>::denorm_min();
    if (Q.col(j).norm() <= eps * A.col(j).norm()) {
      std::cerr << "Gram-Schmidt failed because "
                << "A has (almost) linearly dependent "
                << "columns." << std::endl;
      break;
    } else {
      Q.col(j).normalize();
    }
  }
  // END

  return Q;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
bool testGramSchmidt(unsigned int n) {
  // Orthonormality test
  bool orthonormal = false;
  constexpr double eps = 1e-9;
  Eigen::MatrixXd A, Q;

  // TODO: (1-2.c) Create A, use gram\_schmidt() to compute an
  // orthonormalization of A, call it Q, and return whether it is orthonormal.
  // START
  Eigen::VectorXd linspace = Eigen::VectorXd::LinSpaced(n, 1., n);
  A = Eigen::MatrixXd::Zero(n, n);
  A.rowwise() += 2. * linspace.transpose();
  A.colwise() += linspace;
  Q = gram_schmidt(A);
  const double err =
      (Q.transpose() * Q - Eigen::MatrixXd::Identity(n, n)).norm();
  orthonormal = err < eps;
  // END

  return orthonormal;
}
/* SAM_LISTING_END_1 */

#endif