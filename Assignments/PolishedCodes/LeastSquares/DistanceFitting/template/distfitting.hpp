#ifndef DISTFITTING_HPP
#define DISTFITTING_HPP

#include <Eigen/Sparse>

#include "totriplets.hpp"

/**
 * @brief Initializes the system matrix A
 *
 * @param n number of points
 * @return Eigen::SparseMatrix<double> A
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::SparseMatrix<double> initA(unsigned int n) {
  const unsigned int rows = n * (n - 1) / 2;
  const unsigned int cols = n - 1;
  Eigen::SparseMatrix<double> A(rows, cols);
  // TODO: (3-10.b) Initialize the sparse coefficient matrix for
  // the distance fitting problem
  // START

  // END
  A.makeCompressed();
  return A;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Solves the extended normal equations for the given problem.
 *
 * @param D n times n matrix whose strict upper triangular part contains the
 * distances
 * @return Eigen::VectorXd solution vector
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd solveExtendedNormalEquations(const Eigen::MatrixXd& D) {
  Eigen::VectorXd x;
  const unsigned int n = D.cols();
  const unsigned int m = n * (n - 1) / 2;
  // TODO: (3-10.c) Solve the extended normal equations with a sparse solver
  // that Eigen provides
  // START

  // END
  return x;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solves the normal equations for the given problem.
 *
 * @param D n times n matrix whose strict upper triangular part contains the
 * distances
 * @return Eigen::VectorXd solution vector
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solveNormalEquations(const Eigen::MatrixXd& D) {
  Eigen::VectorXd x;
  const unsigned int n = D.cols();
  const unsigned int m = n * (n - 1) / 2;
  // TODO: (3-10.e) Solve the normal equations by exploiting the
  // Sherman-Morrison-Woodbury formula, i.e. by using your results from the
  // previous subproblem
  // START

  // END
  return x;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Tests the normal equation methods.
 *
 * @param D n times n matrix whose strict upper triangular part contains the
 * distances
 * @return true if the results of both methods agree up to tol
 * @return false otherwise
 */
/* SAM_LISTING_BEGIN_3 */
bool testNormalEquations(const Eigen::MatrixXd& D) {
  constexpr double tol = 1e-9;  // tolerance to check for
  // TODO: (3-10.f) Call your implementations of solveExtendedNormalEquations()
  // and solveNormalEquations() and return true, if their results agree.
  // START

  // END
  return false;
}
/* SAM_LISTING_END_3 */

#endif