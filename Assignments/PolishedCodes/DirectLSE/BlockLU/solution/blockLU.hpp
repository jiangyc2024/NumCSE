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

  // TODO: (3-3.c) (optional) Solve Ry=c using backward substitution. May help
  // when implementing solve_LSE.
  // START
  // Since R is upper triangular, we can solve by backwards substitution
  for (int i = n - 1; i >= 0; --i) {
    y(i) = c(i);
    for (int j = n - 1; j > i; --j) {
      y(i) -= R(i, j) * y(j);
    }
    y(i) /= R(i, i);
  }
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

  // TODO: (3-3.c) Solve the system Ax=b by LU-Decomposition.
  // START
  // Solve Ly = b through forward substitution.
  // Due to the special structure of our L,
  // the first n entries of y are easy:
  y.head(n) = b.head(n);
  // The last element of y is given by $y_n = b_n -
  // u^T\mathbf{R}^{-1}y_{0...n-1}$
  y(n) = b(n) - u.transpose() * solve_R(R, y.head(n));

  // Solve Ux = y by backward substitution.
  // First we build U
  Eigen::MatrixXd U(n + 1, n + 1);
  U << R, v, Eigen::VectorXd::Zero(n).transpose(),
      -u.transpose() * solve_R(R, v);

  // Then we solve Ux = y
  x = solve_R(U, y);
  // \iffalse Latex comment-delimiter so this part doesn't appear
  // in the solution and screws up formatting
  // Note this could be done using less memory by not constructing
  // U, doing the first step of the back substitution "by hand" and then
  // calling solve_R()
  // Eigen::VectorXd x(n+1);
  // x(n) = y(n) / (-u.transpose()*solve_R(R,v));
  // x.head(n) = solve_R(R,y.head(n) - v*x(n));
  //\fi
  // END
  return x;
}
/* SAM_LISTING_END_1 */

#endif