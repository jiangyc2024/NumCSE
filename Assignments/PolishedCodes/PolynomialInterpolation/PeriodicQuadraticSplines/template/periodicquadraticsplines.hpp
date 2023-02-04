#ifndef PERIODICQUADRATICSPLINESHPP
#define PERIODICQUADRATICSPLINESHPP

/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

/* SAM_LISTING_BEGIN_1 */
class ClosedQuadraticSplineCurve {
 public:
  // Constructor
  explicit ClosedQuadraticSplineCurve(const Eigen::Matrix2Xd &p);
  ~ClosedQuadraticSplineCurve() = default;
  ClosedQuadraticSplineCurve() = delete;
  ClosedQuadraticSplineCurve(const ClosedQuadraticSplineCurve &) = delete;
  ClosedQuadraticSplineCurve(const ClosedQuadraticSplineCurve &&) = delete;
  ClosedQuadraticSplineCurve &operator=(const ClosedQuadraticSplineCurve &) =
      delete;
  // Point evaluation operator for sorted parameter arguments
  [[nodiscard]] Eigen::Matrix2Xd curve_points(const Eigen::VectorXd &v) const;
  // Curvature evaluation operator for sorted parameter arguments
  [[nodiscard]] Eigen::VectorXd local_curvatures(
      const Eigen::VectorXd &v) const;
  // Approximate computation of the length of the curve
  [[nodiscard]] double length(double rtol = 1E-6) const;

 private:
  [[nodiscard]] bool checkC1(double tol = 1.0E-4) const;
  // Number of points to be interpolated
  unsigned int n_;
  // Coordinates of interpolated points, duplicate: $\cob{\Vp^0=\Vp^n}$
  Eigen::Matrix2Xd p_;
  // Knot sequence
  Eigen::VectorXd t_;
  // Coefficients in local representation
  Eigen::Matrix2Xd x_;
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
ClosedQuadraticSplineCurve::ClosedQuadraticSplineCurve(
    const Eigen::Matrix2Xd &p)
    : n_(p.cols()), p_(2, p.cols() + 1), t_(p.cols() + 1), x_(2, p.cols()) {
  assert((n_ > 2) && "At least three points have to be supplied");
  assert((n_ % 2 == 1) &&
         "Number of points must be odd!");  // \Label[line]{pqs:1}
  // Save point coordinates, duplicate endpoints
  p_.col(0) = p.col(n_ - 1);
  p_.rightCols(n_) = p;
  // TODO: (5-15.c) Write the constructor.
  // START

  // END
}
/* SAM_LISTING_END_2 */

bool ClosedQuadraticSplineCurve::checkC1(double tol) const {
  bool ok = true;
  // Lengths of knot intervals
  Eigen::VectorXd h = t_.tail(n_) - t_.head(n_);
  double diff = 1.0E-8;
  // Last point
  double tl = 1.0 - diff / h[n_ - 1];
  Eigen::Vector2d left_val = (1.0 - tl) * p_.col(n_ - 1) +
                             h[n_ - 1] * x_.col(n_ - 1) * tl * (1.0 - tl) +
                             tl * p_.col(n_);
  double tr = diff / h[0];
  Eigen::Vector2d right_val = (1.0 - tr) * p_.col(0) +
                              h[0] * x_.col(0) * tr * (1.0 - tr) +
                              tr * p_.col(1);
  Eigen::Vector2d left_slope = (left_val - p_.col(0)) / diff;
  Eigen::Vector2d right_slope = (p_.col(n_) - right_val) / diff;
  std::cout << "CHECK: at p[0]=p[n]: slope diff = "
            << (left_slope - right_slope).norm() << std::endl;
  // First point
  tl = 1.0 - diff / h[0];
  tr = diff / h[0];
  left_val = (1.0 - tl) * p_.col(0) + h[0] * x_.col(0) * tl * (1.0 - tl) +
             tl * p_.col(1);
  right_val = (1.0 - tr) * p_.col(1) + h[1] * x_.col(1) * tr * (1.0 - tr) +
              tr * p_.col(2);
  left_slope = (left_val - p_.col(1)) / diff;
  right_slope = (p_.col(1) - right_val) / diff;
  std::cout << "CHECK: at p[0]=p[n]: slope diff = "
            << (left_slope - right_slope).norm() << std::endl;
  // Internal points
  for (unsigned int k = 1; k < n_ - 1; ++k) {
    tl = 1.0 - diff / h[k];
    tr = diff / h[k + 1];
    left_val = (1.0 - tl) * p_.col(k) + h[k] * x_.col(k) * tl * (1.0 - tl) +
               tl * p_.col(k + 1);
    right_val = (1.0 - tr) * p_.col(k + 1) +
                h[k + 1] * x_.col(k + 1) * tr * (1.0 - tr) + tr * p_.col(k + 2);
    left_slope = (left_val - p_.col(k + 1)) / diff;
    right_slope = (p_.col(k + 1) - right_val) / diff;
    std::cout << "CHECK: at p[0]=p[n]: slope diff = "
              << (left_slope - right_slope).norm() << std::endl;
  }
  return ok;
}

// Computation of points on the curve for many parameter values
// Those are assumed to be sorted
/* SAM_LISTING_BEGIN_3 */
Eigen::Matrix2Xd ClosedQuadraticSplineCurve::curve_points(
    const Eigen::VectorXd &v) const {
  unsigned int N = v.size();
  assert(N > 0);
  // Matrix containing points to be computed
  Eigen::Matrix2Xd s(2, N);
  // Lengths of knot intervals
  const Eigen::VectorXd h = t_.tail(n_) - t_.head(n_);
  // TODO: (5-15.d) Compute the points on the curve for the sorted parameter
  // values v.
  // START

  // END
  return s;
}
/* SAM_LISTING_END_3 */

// Computation of local curvatures of the curve for many parameter values
// Those are assumed to be sorted
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd ClosedQuadraticSplineCurve::local_curvatures(
    const Eigen::VectorXd &v) const {
  unsigned int N = v.size();
  assert(N > 0);
  // Vector for storing the computed curvatures
  Eigen::VectorXd c(N);
  // Lengths of knot intervals
  const Eigen::VectorXd h = t_.tail(n_) - t_.head(n_);
  // TODO: (5-15.e) Compute the local curvature at the parameter values v.
  // START

  // END
  return c;
}
/* SAM_LISTING_END_4 */

// Adaptive computation of the length of the curve by means
// of polygonal approximation up to a given relative tolerance
double ClosedQuadraticSplineCurve::length(double rtol) const {
  double length = 0.0;

  // TODO: (5-15.f) Return an approximation of the length of the curve.
  // START

  // END

  return length;
}

#endif
