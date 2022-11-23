#include <Eigen/Dense>

/* Levinson algorith for solving a Yule-Walker linear system of equations */

/* SAM_LISTING_BEGIN_0 */
void levinson(const Eigen::VectorXd &u, const Eigen::VectorXd &b,
              Eigen::VectorXd &x, Eigen::VectorXd &y) {
  int k = u.size() - 1; // Matrix size - 1
  // Trivial case of $1\times1$ linear sysrtem 
  if (k == 0) {
    x.resize(1); x(0) = b(0);
    y.resize(1); y(0) = u(0);
    return;
  }
  // Vectors holding result of recursive call
  Eigen::VectorXd xk, yk;
  // Recursive call for computing $\cob{\Vx^k}$ and $\cob{\Vy^k}$
  levinson(u.head(k), b.head(k), xk, yk);
  // Coefficient $\cob{\sigma_k}$ from \eqref{eq:ywalg}
  const double sigma = 1 - u.head(k).dot(yk);
  // Update of $\cob{\Vx}$ according to \eqref{eq:ywalg}
  const double t = (b(k) - u.head(k).reverse().dot(xk)) / sigma;
  x = xk - t * yk.head(k).reverse();
  x.conservativeResize(x.size() + 1);
  x(x.size() - 1) = t;
  // Update of vectors $\cob{\Vy^k}$
  double s = (u(k) - u.head(k).reverse().dot(yk)) / sigma;
  y = yk - s * yk.head(k).reverse();
  y.conservativeResize(y.size() + 1);
  y(y.size() - 1) = s;
}
/* SAM_LISTING_END_0 */
