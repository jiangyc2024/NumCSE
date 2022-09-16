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

#endif
