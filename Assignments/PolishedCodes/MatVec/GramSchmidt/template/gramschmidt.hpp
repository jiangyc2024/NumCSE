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

  // END

  return orthonormal;
}
/* SAM_LISTING_END_1 */

#endif