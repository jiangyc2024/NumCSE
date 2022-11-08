#ifndef QRITERATION_HPP
#define QRITERATION_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

Eigen::Vector2d givens(Eigen::Vector2d a) {
  if (a[1] != 0.0) {
    double s;
    double c;
    if (std::abs(a[1]) > std::abs(a[0])) {  // Avoid cancellation/overflow
      const double t = -a[0] / a[1];
      s = 1.0 / std::sqrt(1.0 + t * t);
      c = s * t;
    } else {
      const double t = -a[1] / a[0];
      c = 1.0 / std::sqrt(1 + t * t);
      s = c * t;
    }
    return Eigen::Vector2d(c, s);
  }
  return Eigen::Vector2d(1.0, 0.0);
}

/**
 * @brief Computes the defining vectors of the symmetric, tridiagonal matrix
 * $\cob{\VT'} = \cob{\VQ^\top \VT \VQ}$
 *
 * @param dT diagonal of tridiagonal, symmetric matrix $\cob{\VT}$, size $n$
 * @param uT subdiagonal of tridiagonal, symmetric matrix $\cob{\VT}$, size $n -
 * 1$
 * @return std::pair<Eigen::VectorXd, Eigen::VectorXd> tuple of defining vectors
 * of $\cob{\VT'}$
 */
/* SAM_LISTING_BEGIN_7 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> qrStep(const Eigen::VectorXd &dT,
                                                   const Eigen::VectorXd &uT) {
  const long n = dT.size();
  assert(uT.size() == n - 1);
  // Defining vectors for T' (Written as $T_p$)
  Eigen::VectorXd dT_p = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd uT_p = Eigen::VectorXd::Zero(n - 1);

  // TODO: (3-13.c) Compute the defining vectors d(T') and u(T') of T':= RQ with
  // O(n) asymptotic complexity, given the defining vectors d(T) and U(T) of T =
  // QR
  // START

  // END

  return {dT_p, uT_p};
}
/* SAM_LISTING_END_7 */

#endif