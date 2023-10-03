///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

namespace arrowsys {


using Eigen::VectorXd;

//NOLINTBEGIN(bugprone-easily-swappable-parameters)
inline
/* SAM_LISTING_BEGIN_0 */
VectorXd arrowsys_fast(const VectorXd &d,
                       const VectorXd &c,
                       const VectorXd &b, double alpha,
                       const VectorXd &y) {
  const Eigen::Index n = d.size();
  const VectorXd z = c.array() / d.array();  // \Blue{$\Vz = \VD^{-1}\Vc$}
  const VectorXd w =
      y.head(n).array() / d.array();    // \Blue{$\Vw = \VD^{-1}\Vy_1$}
  const double den = alpha - b.dot(z);  // denominator in \eqref{eqmv:am3}
  // Check for (relatively!) small denominator
  if (std::abs(den) < std::numeric_limits<double>::epsilon() *
                          (std::abs(b.dot(z)) + std::abs(alpha))) {
    throw std::runtime_error("Nearly singular system");
  }
  const double xi = (y(n) - b.dot(w)) / den;
  return (Eigen::VectorXd(n + 1) << w - xi * z, xi).finished();
}
/* SAM_LISTING_END_0 */
//NOLINTEND(bugprone-easily-swappable-parameters)


} // namespace arrowsys