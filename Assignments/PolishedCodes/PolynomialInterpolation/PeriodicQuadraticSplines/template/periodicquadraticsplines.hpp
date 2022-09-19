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

class ClosedQuadraticSplineCurve {
 public:
  // Constructor
  explicit ClosedQuadraticSplineCurve(
      const Eigen::Matrix<double, 2, Eigen::Dynamic> &p);
  ~ClosedQuadraticSplineCurve() = default;
  ClosedQuadraticSplineCurve() = delete;
  ClosedQuadraticSplineCurve(const ClosedQuadraticSplineCurve &) = delete;
  ClosedQuadraticSplineCurve(const ClosedQuadraticSplineCurve &&) = delete;
  ClosedQuadraticSplineCurve &operator=(const ClosedQuadraticSplineCurve &) =
      delete;
  // Point evaluation operator for sorted parameter arguments
  [[nodiscard]] Eigen::Matrix<double, 2, Eigen::Dynamic> curve_points(
      const Eigen::VectorXd &v) const;
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
  Eigen::Matrix<double, 2, Eigen::Dynamic> p_;
  // Knot sequence
  Eigen::VectorXd t_;
  // Coefficients in local representation
  Eigen::Matrix<double, 2, Eigen::Dynamic> x_;
};

ClosedQuadraticSplineCurve::ClosedQuadraticSplineCurve(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &p)
    : n_(p.cols()), p_(2, p.cols() + 1), t_(p.cols() + 1), x_(2, p.cols()) {
  assert((n_ > 2) && "At least three points have to be supplied");
  assert((n_ % 2 == 1) && "Number of points must be odd!"); // \Label[line]{pqs:1}
  // Save point coordinates, duplicate endpoints
  p_.col(0) = p.col(n_ - 1);
  p_.rightCols(n_) = p;
  // START Student code

  // END Student code
}

// Computation of points on the curve for many parameter values
// Those are assumed to be sorted
Eigen::Matrix<double, 2, Eigen::Dynamic>
ClosedQuadraticSplineCurve::curve_points(const Eigen::VectorXd &v) const {
  unsigned int N = v.size();
  assert(N > 0);
  // Matrix containing points to be computed
  Eigen::Matrix<double, 2, Eigen::Dynamic> s(2, N);
  // TO DO: Implement an efficient function to evaluate the spline at sorted
  //        input points v
  // START Student code

  // END Student code
  return s;
}

Eigen::VectorXd ClosedQuadraticSplineCurve::local_curvatures(
    const Eigen::VectorXd &v) const {
  unsigned int N = v.size();
  assert(N > 0);
  // Vector for storing the computed curvatures
  Eigen::VectorXd c(N);
  // Lengths of knot intervals
  const Eigen::VectorXd h = t_.tail(n_) - t_.head(n_);
  // START Student code

  // END Student code
  return c;
}

// Adaptive computation of the length of the curve by means
// of polygonal approximation up to a given relative tolerance
double ClosedQuadraticSplineCurve::length(double rtol) const {
  double length = 0.0;

  return length;
}


#endif
