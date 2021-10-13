
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

std::pair<Eigen::VectorXd, Eigen::VectorXd> qrStep(const Eigen::VectorXd &dT,
                                                   const Eigen::VectorXd &uT) {
  const long n = dT.size();
  assert(uT.size() == n - 1);
  // Defining vectors for T' (Written as T_p)
  Eigen::VectorXd dT_p(n);
  Eigen::VectorXd uT_p(n-1);
  
  // TO DO : Compute the defining vectors d(T') and u(T') of T':= RQ with O(n)  
  // asymptotic complexity, given the defining vectors d(T) and U(T) of T = QR
  // START
  
  // END
  
  return {dT_p, uT_p};
}
