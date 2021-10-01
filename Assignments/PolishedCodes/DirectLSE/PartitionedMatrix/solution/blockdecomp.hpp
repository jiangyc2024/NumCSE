#ifndef BLOCKDECOMP_HPP
#define BLOCKDECOMP_HPP

#include <Eigen/Dense>
#include <cassert>

/**
 * @brief Use efficient implementation A*x = bb
 *
 * @param R MatrixXd is nxn and upper triangular
 * @param v VectorXd is nx1
 * @param u VectorXd is nx1
 * @param bb vector is (n+1)x1 and is stacked (b, \beta)^T =: b
 * @param x solution A*bb = x
 */
/* SAM_LISTING_BEGIN_1 */
void solvelse(const Eigen::MatrixXd& R, const Eigen::VectorXd& v,
              const Eigen::VectorXd& u, const Eigen::VectorXd& bb,
              Eigen::VectorXd& x) {
  // size of R, which is size of u, v, and size of bb is n+1
  const unsigned int n = R.rows();

  // TODO: (3-9.d)
  // i) Use assert() to check that R, v, u, bb all have the appropriate sizes.
  // ii) Use (3-9.b) to solve the LSE.
  // Hint: Use R.triangularView<Eigen::Upper>() to make use of the triangular
  // structure of R.
  // START
  assert(n == R.cols() && n == u.size() && n == v.size() &&
         n + 1 == bb.size() && "Size mismatch!");

  const Eigen::TriangularView<const Eigen::MatrixXd, Eigen::Upper>& triR =
      R.triangularView<Eigen::Upper>();

  // $s$ is the Schur complement and, in this case, is a scalar
  // also $b_s$ is a scalar
  // $snv = s^{-1}$, $b_s$ as in lecture notes
  // $sinvbs = s^{-1}*b_s$
  const double sinv = -1. / u.dot(triR.solve(v));
  const double bs = (bb(n) - u.dot(triR.solve(bb.head(n))));
  const double sinvbs = sinv * bs;

  // Stack the vector $(z, \xi)^T =: x$
  x = Eigen::VectorXd::Zero(n + 1);
  x << triR.solve(bb.head(n) - v * sinvbs), sinvbs;
  // END
}
/* SAM_LISTING_END_1 */

/**
 * @brief Use Eigen's LU-solver to solve Ax = y and check against solvelse()
 *
 * @param R MatrixXd is nxn and upper triangular
 * @param v VectorXd is nx1
 * @param u VectorXd is nx1
 * @param b vector is (n+1)x1 and is stacked $(b, \beta)^T =: b$
 * @param x solution A*bb = x
 * @return true if the result of Eigen's solver conforms with solvelse()
 * @return false otherwise
 */
/* SAM_LISTING_BEGIN_2 */
bool testSolveLSE(const Eigen::MatrixXd& R, const Eigen::VectorXd& v,
                  const Eigen::VectorXd& u, const Eigen::VectorXd& b,
                  Eigen::VectorXd& x) {
  bool areTheSame = false;

  // TODO: (3-9.e)
  // i) Create the system matrix A and solve the LSE, using
  // an Eigen LU-solver. Store the solution in x.
  // ii) Solve the LSE with solvelse(), and calculate the
  // difference between this solution and x. Return true
  // if and only if the norm of this difference is close
  // enough to zero.
  // START
  // size of R, which is size of u, v, and size of bb is n+1
  const unsigned int n = R.rows();
  assert(n == R.cols() && n == u.size() && n == v.size() && n + 1 == b.size() &&
         "Size mismatch!");

  Eigen::MatrixXd A(n + 1, n + 1);
  // Build matrix using blocks
  A << R, v, u.transpose(), 0;

  // Solve system using LU decpomposition
  x = A.partialPivLu().solve(b);

  // Solve system using solvelse()
  Eigen::VectorXd y;
  solvelse(R, v, u, b, y);

  // We choose an acceptable tolerance.
  constexpr double tolerance = 1e-8;
  areTheSame = ((x - y).norm() < tolerance * b.norm());

  // END
  return areTheSame;
}
/* SAM_LISTING_END_2 */

#endif