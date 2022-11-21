/**
 * @file contourplot_test.cc
 * @brief NPDE homework ContourPlot code
 * @author Oliver Rietmann
 * @date 25.03.2021
 * @copyright Developed at ETH Zurich
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"
#include <Eigen/Core>

struct TestData {
  TestData() {
    T = 6;
    y0 << 2, 0;
  }

  double T;
  Eigen::Vector2d y0;
};

TestData data;

TEST_SUITE("ContourPlot") {
  TEST_CASE("Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePoints" *
            doctest::description("bl")) {
    auto gradF = [](Eigen::Vector2d x) -> Eigen::Vector2d {
      return Eigen::Vector2d(x(0), 2.0 * x(1));
    };
    auto sol = computeIsolinePoints(gradF, data.y0, data.T);
    auto stud = computeIsolinePoints_TEST(gradF, data.y0, data.T);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEgg" *
            doctest::description("abc")) {
    auto sol = crookedEgg();
    auto stud = crookedEgg_TEST();
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePointsDQ" *
            doctest::description("abc")) {
    auto sol = computeIsolinePointsDQ(gradF, data.y0, data.T);
    auto stud = computeIsolinePointsDQ_TEST(gradF, data.y0, data.T);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }
}
