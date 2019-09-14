///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
VectorXd arrowsys_fast(const VectorXd &d, const VectorXd &c, const VectorXd &b,
                       const double alpha, const VectorXd &y) {
  int n = d.size();
  VectorXd z = c.array() / d.array();         // \Blue{$\Vz = \VD^{-1}\Vc$}
  VectorXd w = y.head(n).array() / d.array(); // \Blue{$\Vw = \VD^{-1}\Vy_1$}
  double xi = (y(n) - b.dot(w)) / (alpha - b.dot(z));
  return (VectorXd(n + 1) << w - xi * z, xi).finished();
}
/* SAM_LISTING_END_0 */
 
