#ifndef TRIDIAGLEASTSQUARES_HPP
#define TRIDIAGLEASTSQUARES_HPP

#include <Eigen/Dense>
#include <cassert>

/**
 * @brief Minimize |Tz - c|_2 wrt (alpha,beta), where T is triagonal,
 * T.diagonal()=alpha, T.diagonal(+/-1)=beta.
 *
 * @param z An $n$-dimensional vector of data for lhs of Tz=c.
 * @param c An $n$-dimensional vector of data for rhs of Tz=c.
 * @return Eigen::VectorXd The parameters of T, x=(alpha,beta), that solve the
 * least squares problem Tz=c.
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd lsqEst(const Eigen::VectorXd& z, const Eigen::VectorXd& c) {
  assert(z.size() == c.size() && "z and c must have same size");
  Eigen::VectorXd x = Eigen::VectorXd::Zero(2);  // $(\alpha,\beta)$.

  // TODO: (3-2.c) Solve eq. (3.2.1) using the method of normal equations.
  // START

  // Reformulating $T_{\alpha,\beta}z=c$ as $Ax=b$ with $x=(\alpha,\beta)^T$,
  // we have $b=c$, and $A$ defined as follows.
  const unsigned int n = z.size();
  Eigen::MatrixXd A(n, 2);
  // Initialize all entries of A:
  A.col(0) = z;
  A(0, 1) = 0;
  A.col(1).tail(n - 1) = z.head(n - 1);
  // Add remaining terms:
  A.col(1).head(n - 1) += z.tail(n - 1);
  // Normal equation linear system
  Eigen::MatrixXd lhs = A.transpose() * A;  // Left-hand side
  Eigen::VectorXd rhs = A.transpose() * c;  // Right-hand side
  // Least squares estimate of $(\alpha,\beta)$ as solution
  // of normal equations. 
  x = lhs.lu().solve(rhs);
  // END
  return x;
}
/* SAM_LISTING_END_0 */

#endif
