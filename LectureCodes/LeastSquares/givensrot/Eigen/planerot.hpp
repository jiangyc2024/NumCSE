///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>
#include <cmath>


namespace planerot {


using Eigen::Vector2d;
using Eigen::Matrix2d;

inline
/* SAM_LISTING_BEGIN_0 */
// plane (2D) Givens rotation avoiding cancellation
// Computes \com{orthogonal} $\cob{\VG\in\bbR^{2,2}}$ with $\cob{\VG^\top\Va = \twovec{r}{0} =: \Vx}$, $\cob{r=\pm \N{\Va}_2}$
void planerot(const Eigen::Vector2d& a, Eigen::Matrix2d& G,
              Eigen::Vector2d& x) {
  int sign{1};
  if (a[1] != 0.0) {
    double s = NAN; // s $\leftrightarrow$ $\cob{\sigma}$
    double c = NAN; // c $\leftrightarrow$ $\cob{\gamma}$
    if (std::abs(a[1]) > std::abs(a[0])) {  // Avoid overflow
      const double t = -a[0] / a[1];
      s = 1.0 / std::sqrt(1.0 + t * t);
      c = s * t;
      sign = -1;
    } else {
      const double t = -a[1] / a[0];
      c = 1.0 / std::sqrt(1.0 + t * t);
      s = c * t;
    }
    G << c, s, -s, c;  // Form $2\times 2$ Givens rotation matrix
  } else {
    G.setIdentity();
  }
  x << (sign * a.norm()), 0.0;
}
/* SAM_LISTING_END_0 */

inline
/* SAM_LISTING_BEGIN_1 */
// plane (2D) Givens rotation
// unstable implementation !
void planerot_inst(const Vector2d& a, Matrix2d& G, Vector2d& x) {
  if (a(1) != 0) {
    const double r = a.norm();
    G.row(0) = a.transpose() / r;
    G(1, 0) = -a(1) / r;
    G(1, 1) = a(0) / r;
    x(0) = r;
    x(1) = 0;
  } else {
    G.setIdentity();
  }
}
/* SAM_LISTING_END_1 */

inline
/* SAM_LISTING_BEGIN_2 */
// Computing entries \Blue{$\gamma$} and \Blue{$\sigma$}
// real of Givens matrix \eqref{eq:givens}
std::pair<double, double> planerot(const Vector2d& a) {
  double s = 0;
  double c = 1;
  // Stable computation avoiding cancellation
  if (a(1) != 0.0) {
    if (std::abs(a(1)) > std::abs(a(0))) {
      const double t = -a(0) / a(1);
      s = 1.0 / std::sqrt(1.0 + t * t);
      c = s * t;
    } else {
      const double t = -a(1) / a(0);
      c = 1.0 / std::sqrt(1.0 + t * t);
      s = c * t;
    }
  }
  return {c, s};
}
/* SAM_LISTING_END_2 */


} // namespace planerot