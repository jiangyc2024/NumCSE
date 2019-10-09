// TO DO: Implement linReg() and expFit().
// Note: Your code will not compile until after you have defined both functions.
// START

#include <Eigen/Dense>

using namespace Eigen;

/* @brief Solve the linear regression problem (fitting a line to data)
 * for data points $(t_i,y_i)$, $i = 1,\dots,n$
 * @param[in] t An $n$-dimensional vector containing one side of input data
 * @param[in] y An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(x_1,x_2)$, intercept and slope of the line fitted
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd linReg(const VectorXd &t, const VectorXd &y) {
  assert( t.size() == y.size() && "t and y must have same size" );

  // Coefficient matrix of overdetermined linear system.
  MatrixXd A(t.size(),2);
  A.col(0).setOnes();
  A.col(1) = t;
  
  // Normal equation
  MatrixXd lhs = A.transpose() * A; // Left-hand side
  VectorXd rhs = A.transpose() * y; // Right-hand side
  
  // Intercept and slope.
  VectorXd x = lhs.lu().solve(rhs);
  return x;
}
/* SAM_LISTING_END_0 */

/* @brief Solve the linearized exponential problem
 * for data points $(t_i,y_i)$, $i = 1,\dots,n$
 * @param[in] t An $n$-dimensional vector containing one side of input data
 * @param[in] y An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(x_1,x_2)$, intercept and slope of the line fitted to the linearized problem
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd expFit(const VectorXd &t, const VectorXd &y) {
  assert( t.size() == y.size() && "t and y must have same size");
  
  // Transform to a linear equation
  // b := log(y) = x(0) + x(1)*t
  VectorXd b = y.array().log().matrix();
  // Solve, and transform back to the nonlinear equation
  // y = exp(x(0)) * exp(x(1)*t) = alpha*exp(beta*t)
  VectorXd x = linReg(t, b);
  x(0) = exp(x(0));
  
  return x;
}
/* SAM_LISTING_END_1 */

// END