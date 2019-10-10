// TO DO: Implement lsqEst().
// Note: Your code will not compile until after you have defined the function.
// Note: The test anticipates the output in the order (alpha, beta), and not (beta, alpha).
// START
#include <Eigen/Dense>
using namespace Eigen;

/* @brief Minimize  |Tz - c|_2 wrt (alpha,beta), where T is triagonal, T.diagonal()=alpha, T.diagonal(+/-1)=beta.
 * @param[in] z An $n$-dimensional vector of data for lhs of Tz=c.
 * @param[in] c An $n$-dimensional vector of data for rhs of Tz=c.
 * @param[out] x The parameters of T, x=(alpha,beta), that solve the least squares problem Tz=c.
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd lsqEst(const VectorXd &z, const VectorXd &c) {
  assert( z.size() == c.size() && "z and c must have same size" );
  
  // Reformulating $Tz=c$ as $Ax=b$ with $x=(\alpha,\beta)^T$,
  // we have $b=c$, and $A$ defined as follows.
  int n = z.size();
  MatrixXd A(n,2);
  // Initialize all entries of A:
  A.col(0) = z;
  A(0,1) = 0;
  A.col(1).tail(n-1) = z.head(n-1);
  // Add remaining terms:
  A.col(1).head(n-1) += z.tail(n-1);
  
  // Normal equation
  MatrixXd lhs = A.transpose() * A; // Left-hand side
  VectorXd rhs = A.transpose() * c; // Right-hand side

  // Intercept and slope.
  VectorXd x = lhs.lu().solve(rhs);
  return x;
}
/* SAM_LISTING_END_0 */

//END