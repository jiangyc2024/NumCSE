///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd arrowsys_fast(const Eigen::VectorXd &d,
                              const Eigen::VectorXd &c,
                              const Eigen::VectorXd &b, double alpha,
                              const Eigen::VectorXd &y) {
  int n = d.size();
  Eigen::VectorXd z = c.array() / d.array();  // \Blue{$\Vz = \VD^{-1}\Vc$}
  Eigen::VectorXd w =
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
