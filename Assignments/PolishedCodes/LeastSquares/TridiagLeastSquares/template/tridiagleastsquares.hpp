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
  Eigen::VectorXd x = Eigen::VectorXd::Zero(2);  // $(alpha,\beta)$.

  // TODO: (3-2.c) Solve eq. (3.2.1) using the method of normal equations.
  // START

  // END
  return x;
}
/* SAM_LISTING_END_0 */

#endif