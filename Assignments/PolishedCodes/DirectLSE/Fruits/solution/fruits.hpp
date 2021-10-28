#ifndef FRUITS_HPP
#define FRUITS_HPP

#include <Eigen/Dense>
#include <cassert>

/**
 * @brief Solve the LSE arising from
 * the problem description to determine the
 * prices of the fruits using Eigen
 *
 * @return Eigen::VectorXd the fruit prices
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd fruitPrice() {
  Eigen::VectorXd x = Eigen::VectorXd::Zero(6);

  // TODO: (2-1.b) Determine the fruit prices using Eigen.
  // START
  // Initialize the matrix
  Eigen::MatrixXd A(6, 6);
  // See "Advanced Initialization" in the Eigen doc
  A << 3, 1, 7, 2, 0, 0, 2, 0, 0, 0, 0, 1, 1, 0, 0, 0, 3, 0, 0, 5, 0, 1, 0, 0,
      0, 0, 0, 2, 0, 1, 0, 1, 20, 0, 0, 0;

  // Initialize the vector
  Eigen::VectorXd b(6);
  b << 11.10, 17.00, 6.10, 5.25, 12.50, 7.00;

  // Initialize the solver; other decompositions also possible
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver;
  // Decompose A
  solver.compute(A);
  // Solve the LSE
  x = solver.solve(b);

  // Notice how we split up initialization, decompositon and solving.
  // This means we can use the same solver for multiple matrices and the
  // same decompositions for multiple LSE with the same matrix.
  // These three steps can be combined:
  Eigen::VectorXd x2 = A.colPivHouseholderQr().solve(b);

  // The result is the same:
  assert((x - x2).norm() < 1e-10);
  // END
  return x;
}
/* SAM_LISTING_END_1 */

#endif
