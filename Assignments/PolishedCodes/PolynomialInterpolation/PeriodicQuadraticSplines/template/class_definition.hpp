/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

 #ifndef CLASSDEFHPP
 #define CLASSDEFHPP

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
  Eigen::Matrix<double, 2, Eigen::Dynamic> curve_points(
      const Eigen::VectorXd &v) const;
 private:
  // Number of points to be interpolated
  unsigned int n_;
  // Coordinates of interpolated points, duplicate: $\cob{\Vp^0=\Vp^n}$
  Eigen::Matrix<double, 2, Eigen::Dynamic> p_;
  // Knot sequence
  Eigen::VectorXd t_;
  // Coefficients in local representation
  Eigen::Matrix<double, 2, Eigen::Dynamic> x_;
};

#endif
