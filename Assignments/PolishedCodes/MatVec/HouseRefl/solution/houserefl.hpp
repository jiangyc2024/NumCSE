#ifndef HOUSEREFL_HPP
#define HOUSEREFL_HPP

#include <Eigen/Dense>

/* \brief Compute an ONB of the space orthogonal to $v$
 * @param v vector $\in \mathbb{R}^n \setminus \{ 0 \}$
 * @param Z matrix $\in \mathbb{R}^{n-1 \times n}$
 */
/* SAM_LISTING_BEGIN_1 */
void houserefl(const Eigen::VectorXd &v, Eigen::MatrixXd &Z)
{

  // TODO: (1-5.a) Implement the householder refelction algorithm described in
  // the pseudo code 1.5.1 to get the ONB of the space orthogonal to $v$.

  // START
  unsigned int n = v.size();
  Eigen::VectorXd w = v.normalized();
  Eigen::VectorXd u=w;
  u(0) += 1;
  Eigen::VectorXd q = u.normalized();
  Eigen::MatrixXd X = Eigen::MatrixXd::Identity(n, n) - 2*q*q.transpose();
  Z = X.rightCols(n-1);
  // END

}
/* SAM_LISTING_END_1 */

#endif
