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

namespace zerosquadpolstab {


using Eigen::VectorXd;
using Eigen::Vector2d;

inline
/* SAM_LISTING_BEGIN_0 */
//! \cpp function computing the zeros of a quadratic polynomial
//! $\xi\to \xi^2+\alpha\xi+\beta$ by means of the familiar discriminant
//! formula $\xi_{1,2} = \frac{1}{2}(-\alpha\pm\sqrt{\alpha^2-4\beta})$.
//! This is a stable implementation based on Vieta's theorem.
//! The zeros are returned in a column vector
Eigen::VectorXd zerosquadpolstab(double alpha, double beta) {
  Eigen::Vector2d z(2);
  double D = std::pow(alpha, 2) - 4 * beta; // discriminant
  if (D >= 0) {
    double wD = std::sqrt(D);
    // Use discriminant formula only for zero far away from $0$
    // in order to \com{avoid cancellation}. For the other zero
    // use Vieta's formula.
    if (alpha >= 0) {
      double t = 0.5 * (-alpha - wD); // \Label[line]{zqs:11}
      z << t, beta / t;
    } else {
      double t = 0.5 * (-alpha + wD); // \Label[line]{zqs:12}
      z << beta / t, t;
    }
  }
  else {
    throw std::runtime_error( "no real zeros" );
  }
  return z;
}
/* SAM_LISTING_END_0 */


} //namespace zerosquadpolstab
