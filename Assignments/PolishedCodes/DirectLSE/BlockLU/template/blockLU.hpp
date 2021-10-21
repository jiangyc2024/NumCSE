#ifndef BLOCKLU_HPP
#define BLOCKLU_HPP

#include <Eigen/Dense>
#include <cassert>

/**
 * @brief Solve the system Ry=c for the upper triangular matrix R. This could
 * help you in your implementation of solve_LSE().
 *
 * @param R nxn regular, upper triangular matrix
 * @param c n dim vector
 * @return Eigen::VectorXd n dim result vector
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solve_R(const Eigen::MatrixXd& R, const Eigen::VectorXd& c) {
  const unsigned int n = R.rows();
  assert(n == R.cols() && n == c.size() && "Input dimensions must agree");
  // initialize the return value
  Eigen::VectorXd y(n);

  // TODO: (optional) Solve Ry=c using backward substitution. May help
  // when implementing solve\_LSE.
  // START

  // END
  return y;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Solve the System Ax=b
 *        for A << R,              v,
 *                 u.transpose(),  0;
 *
 * @param R nxn regular, upper triangular matrix
 * @param v n dim vector
 * @param u n dim vector
 * @param b n+1 dim vector
 * @return Eigen::VectorXd n+1 dim result vector
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd solve_LSE(const Eigen::MatrixXd& R, const Eigen::VectorXd& v,
                          const Eigen::VectorXd& u, const Eigen::VectorXd& b) {
  const unsigned int n = R.rows();
  assert(R.cols() == n && "R has to be square");
  assert(n == v.size() && n == u.size() && n + 1 == b.size() &&
         "Input dimensions must agree");
  Eigen::VectorXd y(n + 1), x(n + 1);

  // TODO: (2-3.d) Solve the system Ax=b by LU-Decomposition.
  // START

  // END
  return x;
}
/* SAM_LISTING_END_1 */

#endif
