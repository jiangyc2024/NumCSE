/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: 2019
 */

#include <Eigen/Dense>
#include <Eigen/QR>

/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTION, typename JACOBIAN>
Eigen::VectorXd gn(const Eigen::VectorXd &init, FUNCTION &&F, JACOBIAN &&J,
            double rtol = 1.0E-6, double atol = 1.0E-8) {
  Eigen::VectorXd x = init; // Vector for iterates $\cob{\Vx^{(k)}}$
  // Vector for Gauss-Newton correction $\cob{\Vs}$
  Eigen::VectorXd s = J(x).householderQr().solve(F(x)); // \Label[line]{gn:2}
  x = x - s;
  // A posteriori termination based on absolute and relative tolerances
  while ((s.norm() > rtol * x.norm()) && (s.norm() > atol)) { // \Label[line]{gn:term}
    s = J(x).householderQr().solve(F(x)); // \Label[line]{gn:5}
    x = x - s;
  }
  return x;
}
/* SAM_LISTING_END_1 */
