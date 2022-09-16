/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

 #ifndef PERIODICQUADRATICSPLINESHPP
 #define PERIODICQUADRATICSPLINESHPP

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include "class_definition.hpp"

// Computation of points on the curve for many parameter values
// Those are assumed to be sorted
/* SAM_LISTING_BEGIN_3 */
Eigen::Matrix<double, 2, Eigen::Dynamic>
ClosedQuadraticSplineCurve::curve_points(const Eigen::VectorXd &v) const {
  unsigned int N = v.size();
  assert(N > 0);
  // Matrix containing points to be computed
  Eigen::Matrix<double, 2, Eigen::Dynamic> s(2, N);
  // Lengths of knot intervals
  const Eigen::VectorXd h = t_.tail(n_) - t_.head(n_);
  // START Student code
  // Run through all the provided parameter values
  unsigned int knot_idx = 0;
  for (unsigned int k = 0; k < N; ++k) {
    assert(((v[k] >= 0.0) && (v[k] < 1.0)) && "Out of range!");
    if (k > 0) {
      assert(v[k] >= v[k - 1] && "Parameter values not sorted!");
    }
    // Find right knot interval: knot\_idx stores index of knot to the right of
    // current parameter value
    while ((knot_idx < n_) && (v[k] >= t_[knot_idx])) {
      knot_idx++;
    }
    assert(knot_idx > 0);
    const double tau = (v[k] - t_[knot_idx - 1]) / h[knot_idx - 1];
    s.col(k) = ((1.0 - tau) * p_.col(knot_idx - 1)) +
               (h[knot_idx - 1] * tau * (1 - tau)) * x_.col(knot_idx - 1) +
               (tau * p_.col(knot_idx));
  }
  // END Student code
  return s;
}
/* SAM_LISTING_END_3 */
#endif
