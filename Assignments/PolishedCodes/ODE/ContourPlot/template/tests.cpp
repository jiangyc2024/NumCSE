#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Core>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    T = 6;
    y0 << 2, 0;
  }

  double T;
  Eigen::Vector2d y0;
};

TestData data;

Eigen::Vector2d gradF(Eigen::Vector2d x) {
  return Eigen::Vector2d(x(0), 2.0 * x(1));
}

double F(Eigen::Vector2d x) { return 0.5 * x(0) * x(0) + x(1) * x(1); }

TEST_SUITE("ContourPlot") {
  TEST_CASE("Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePoints" *
            doctest::description("Checking all the points")) {
    auto sol = computeIsolinePoints(gradF, data.y0, data.T);
    auto stud = computeIsolinePoints_TEST(gradF, data.y0, data.T);
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEgg" *
            doctest::description("Checking all the points")) {
    auto sol = crookedEgg();
    auto stud = crookedEgg_TEST();
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePointsDQ" *
            doctest::description("Checking all the points")) {
    auto sol = computeIsolinePointsDQ(F, data.y0, data.T);
    auto stud = computeIsolinePointsDQ_TEST(F, data.y0, data.T);
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }
}
