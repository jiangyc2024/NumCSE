#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

// includes for test data
#include "notaknotcubicspline.hpp"
#include <Eigen/Dense>
#include <cmath>

TEST_SUITE("Not-a-knot Cubic Splines") {

  TEST_CASE(
      "NotAKnotCubicSpline()" *
      doctest::description("Testing if eval() works as expected after invoking "
                           "the constructor (Not testing efficiency)")) {
    unsigned N = 4;
    Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(N, -1, 2));
    Eigen::VectorXd y(4);
    y << 0, 1, 4, 15;
    NotAKnotCubicSpline spline(t, y);
    CHECK(std::abs(spline.eval(1.5) - 8.125) ==
          doctest::Approx(0.).epsilon(1e-10));
    CHECK(std::abs(spline.eval(.5) - 1.875) ==
          doctest::Approx(0.).epsilon(1e-10));
  }

  TEST_CASE("evalDerivative()" * doctest::description("Testing the output")) {
    unsigned N = 4;
    Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(N, -1, 2));
    Eigen::VectorXd y(4);
    y << 0, 1, 4, 15;
    NotAKnotCubicSpline spline(t, y);
    CHECK(std::abs(spline.evalDerivative(1.5) - 10.75) ==
          doctest::Approx(0.).epsilon(1e-10));
    CHECK(std::abs(spline.evalDerivative(.5) - 2.75) ==
          doctest::Approx(0.).epsilon(1e-10));
  }
}
