/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#define _USE_MATH_DEFINES

#ifdef NICEBACKTRACE
#include "backtrace.hpp"
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include "periodicquadraticsplines.hpp"

#ifdef NICEBACKTRACE
#include "backtrace.hpp"
#endif

#define SMART_LSE_SOLUTION

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "NumCSE code for close quadratic spline curve" << std::endl;

  {
    // clang-format off
    const Eigen::Matrix<double, 2, Eigen::Dynamic> p =
        (Eigen::Matrix<double, 2, Eigen::Dynamic>(2, 5) <<
	 0, 1, 1, 0.5, 0,
	 0, 0, 1, 1.5, 1)
            .finished();
    // clang-format on
    std::cout << "\nClosedQuadraticSplineCurve object created using points \n \n" << p.transpose() << std::endl;
    ClosedQuadraticSplineCurve s(p);

    // Output points on spline curve
    unsigned int N = 10;
    auto pts{s.curve_points(Eigen::VectorXd::LinSpaced(N, 0.0, N - 1) / N)};
    std::cout << "\nEvatuating Spline at equispaced pts = \n" << std::endl;
    for (unsigned int k = 0; k < pts.cols(); ++k) {
      std::cout << pts.col(k).transpose() << std::endl;
    }
    std::cout << p.col(4).transpose() <<  std::endl;
  }
  return 0;
}
