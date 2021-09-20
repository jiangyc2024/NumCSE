#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    f = [](Eigen::Vector2d x) {
      return std::pow(x(0), 2) + 2 * std::pow(x(1), 4);
    };
    c = 2;
  }
  std::function<double(Eigen::Vector2d)> f;
  double c;
};

TestData data;

TEST_SUITE("Level Set") {
  TEST_CASE("Eigen::Vector2d pointLevelSet" *
            doctest::description("pointLevelSet")) {
    Eigen::Vector2d ptest = {1.28718850581117, 0.643594252905583};
    Eigen::Vector2d d = {2, 1};
    Eigen::Vector2d x0 = {std::sqrt(2), 0};
    Eigen::Vector2d p = pointLevelSet_TEST(data.f, d, data.c, x0, 1e-10, 1e-16);
    CHECK((p - ptest).norm() == doctest::Approx(0).epsilon(1e-14));
  }

  TEST_CASE("double areaLevelSet" * doctest::description("areaLevelSet")) {
    constexpr unsigned int n = 8;
    const double area = areaLevelSet_TEST(data.f, n, data.c);
    constexpr double areatest = 4.26647319712633;
    CHECK(area == doctest::Approx(areatest).epsilon(1e-14));
  }
}
