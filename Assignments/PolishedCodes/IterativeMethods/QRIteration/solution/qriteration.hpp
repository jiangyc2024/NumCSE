#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace Eigen;

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

/* SAM_LISTING_BEGIN_7 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> qrStep(const Eigen::VectorXd &dT,
                                                   const Eigen::VectorXd &uT) {
  const long n = dT.size();
  assert(uT.size() == n - 1);
  // Defining vectors for T' (Written as $T_p$)
  Eigen::VectorXd dT_p(n);
  Eigen::VectorXd uT_p(n-1);
  
  // TO DO : Compute the defining vectors d(T') and u(T') of T':= RQ with O(n)  
  // asymptotic complexity, given the defining vectors d(T) and U(T) of T = QR
  // START
  // Auxiliary 4xn matrix \cob{$\VU\in\bbR^{4,n}$}
  Eigen::MatrixXd U{Eigen::MatrixXd::Zero(4, n)};
  U.block(1, 1, 1, n - 1) = uT.transpose();
  U.block(2, 0, 1, n) = dT.transpose();
  U.block(3, 0, 1, n - 1) = uT.transpose();
  // Stage I: QR decomposition by $\cob{n-1}$ Givens rotations
  std::vector<Eigen::Matrix2d> givensrots; // Store rotations!
  for (long i = 0; i < n - 1; ++i) {
    const Eigen::Vector2d g = givens(Eigen::Vector2d(U(2, i), U(3, i)));
    const Eigen::Matrix2d G{
      (Eigen::Matrix2d(2, 2) << g[0], -g[1], g[1], g[0]).finished()};
    U.block(2, i, 2, 1).applyOnTheLeft(G);
    U.block(1, i + 1, 2, 1).applyOnTheLeft(G);
    if (i < n - 2) {
      U.block(0, i + 2, 2, 1).applyOnTheLeft(G);
    }
    givensrots.push_back(G);
  }
  // At this point $\cob{\VU}$ encodes the R-factor of the QR-decomposition of
  // $\cob{\VT}$. In particular, the bottom row of $\cob{\VU}$ should contain
  // only entries that are "numerically zero''.

  // Stage II: Apply transposed rotations from right to $\cob{\VR}$ according to
  // formulas \prbeqref{eq:II}
  for (long i = 1; i < n; ++i) {
    const Eigen::Matrix2d &G{givensrots[i - 1]};  // An alias
    Eigen::Vector2d tmp;
    tmp = G * Eigen::Vector2d(U(3, i - 1), U(2, i));
    U(3, i - 1) = tmp[0];
    U(2, i) = tmp[1];
    tmp = G * Eigen::Vector2d(U(2, i - 1), U(1, i));
    U(2, i - 1) = tmp[0];
    U(1, i) = tmp[1];
    tmp = G * Eigen::Vector2d(U(1, i - 1), U(0, i));
    U(1, i - 1) = tmp[0];
    U(0, i) = tmp[1];
  }
  // Extract diagonal and sub-/super-diagonal of $\cob{\VT'}$.'
  dT_p = U.block(2, 0, 1, n).transpose();
  uT_p = U.block(1, 1, 1, n - 1).transpose();
  // END
  
  return {dT_p, uT_p};
}
/* SAM_LISTING_END_7 */
