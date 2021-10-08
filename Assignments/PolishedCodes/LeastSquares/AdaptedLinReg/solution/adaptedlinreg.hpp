#include <Eigen/Dense>

/* @brief Solve the linear regression problem (fitting a line to data)
 * for data points $(t_i,y_i)$, $i = 1,\dots,n$
 * @param[in] t An $n$-dimensional vector containing one side of input data
 * @param[in] y An $n$-dimensional vector containing the other side of input
 * data
 * @param[out] x The vector of parameters $(x_1,x_2)$, intercept and slope of
 * the line fitted
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd linReg(const Eigen::VectorXd &t, const Eigen::VectorXd &y) {
  // TODO: (3-1.d) Use the method of normal equations to solve the 1D linear
  // regression problem in least-square sense.
  // Note: The tests anticipate the outputs in the order (alpha, beta), and not
  // (beta, alpha).
  
  // START
  assert(t.size() == y.size() && "t and y must have same size");

  // Coefficient matrix of overdetermined linear system.
  Eigen::MatrixXd A(t.size(), 2);
  A.col(0).setOnes();
  A.col(1) = t;

  // Normal equation
  Eigen::MatrixXd lhs = A.transpose() * A;  // Left-hand side
  Eigen::VectorXd rhs = A.transpose() * y;  // Right-hand side

  // Intercept and slope.
  Eigen::VectorXd x = lhs.lu().solve(rhs);
  
  // END

  return x;
}
/* SAM_LISTING_END_0 */

/* @brief Solve the linearized exponential problem
 * for data points $(t_i,y_i)$, $i = 1,\dots,n$
 * @param[in] t An $n$-dimensional vector containing one side of input data
 * @param[in] y An $n$-dimensional vector containing the other side of input
 * data
 * @param[out] x The vector of parameters $(x_1,x_2)$, intercept and slope of
 * the line fitted to the linearized problem
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd expFit(const Eigen::VectorXd &t, const Eigen::VectorXd &y) {
  // TODO: (3-1.e) Implement least square estimate of alpha and beta using
  // the previously implemented function linReg().
  // Note: You don't need to have implemented linReg() to solve this subproblem.
  // The tests will compile and use the master solution for linReg().
  // Note: The tests anticipate the outputs in the order (alpha, beta), and not
  // (beta, alpha).
  
  // START
  assert(t.size() == y.size() && "t and y must have same size");

  // Transform to a linear equation
  // b := log(y) = x(0) + x(1)*t
  // We use array() for coefficient-wise log() and matrix() to transform it
  // back to a vector.
  Eigen::VectorXd b = y.array().log().matrix();

  // Solve and transform back to the nonlinear equation
  // y = exp(x(0)) * exp(x(1)*t) = alpha*exp(beta*t)
  Eigen::VectorXd x = linReg(t, b);
  x(0) = exp(x(0));  // Eigen's coefficient-wise exp()

  // END

  return x;
}
/* SAM_LISTING_END_1 */

// END
