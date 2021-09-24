////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>

/* \brief Performs Gram-Schidt orthonormalization
 * Given a matrix $\mathbf{A}$ of linearly independent columns,
 * returns the result of a Gram-Schmidt orthonormalization.
 * Unstable GS algorithm: output is prone to cancellation issues.
 * @param $\mathbf{A}$ Matrix of linearly independent columns
 * \return Matrix with ONB of $span(a_1, \cdots, a_n)$ as columns
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd gram_schmidt(const Eigen::MatrixXd &A) {
  // We create a matrix Q with the same size and data of A
  Eigen::MatrixXd Q(A);

  // The first vector just gets normalized
  Q.col(0).normalize();

  // TO DO: (1-2.b) Implement the gram_schmidt procedure by iterating over all
  // other columns of A. 
  
  // START

  // END

  return Q;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
double orthogonality_test() {
  // Orthonormality test
  double err = 1;
  unsigned int n = 9;
  Eigen::MatrixXd A, Q;
  A = Eigen::MatrixXd::Random(n, n);

  // TO DO: (1-2.c) Use gram_schmidt() to compute an orthonormalization of A,
  // call it Q, and let err measure "how far Q is from being orthonormal".
  
  // START
  
  // END

  return err;
}
/* SAM_LISTING_END_1 */
