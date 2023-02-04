#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>
#define PI M_PI

struct TestData {
  TestData() {
    auto c = [](double t) -> Eigen::Vector2d {
      Eigen::Vector2d ct;
      ct << std::cos(2 * PI * t) + 2. / 3. * std::cos(4 * PI * t),
          3. / 2. * std::sin(2 * PI * t);
      return ct;
    };
    n = 50;
    neval = 33;
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(n, 0, 1.0 - 1.0 / n);
    Sigma.resize(n);
    for (unsigned int i = 0; i < n; i++) {
      Sigma[i] = c(t(i));
    }
    x.resize(neval);
    x = Eigen::VectorXd::LinSpaced(neval, 0, 1);
  }
  Eigen::VectorXd x;
  std::vector<Eigen::Vector2d> Sigma;
  std::size_t n, neval;
};

TestData data;

TEST_SUITE("adaptiveintp") {
  TEST_CASE("std::vector<Eigen::Vector2d> closedPolygonalInterpolant" *
            doctest::description("closedPolygonalInterpolant()")) {
    std::vector<Eigen::Vector2d> sol =
        closedPolygonalInterpolant(data.Sigma, data.x);
    std::vector<Eigen::Vector2d> stud =
        closedPolygonalInterpolant_TEST(data.Sigma, data.x);
    REQUIRE(sol.size() == stud.size());
    double err = 0;
    for (unsigned int i = 0; i < sol.size(); i++) {
      err += (sol[i] - stud[i]).norm();
    }
    CHECK(err == doctest::Approx(0.).epsilon(1e-9));
  }

  TEST_CASE("std::vector<Eigen::Vector2d> closedHermiteInterpolant" *
            doctest::description("closedHermiteInterpolant()")) {
    std::vector<Eigen::Vector2d> sol =
        closedHermiteInterpolant(data.Sigma, data.x);
    std::vector<Eigen::Vector2d> stud =
        closedHermiteInterpolant_TEST(data.Sigma, data.x);
    REQUIRE(sol.size() == stud.size());
    double err = 0;
    for (unsigned int i = 0; i < sol.size(); i++) {
      err += (sol[i] - stud[i]).norm();
    }
    CHECK(err == doctest::Approx(0.).epsilon(1e-9));
  }

  // clang-format off
  TEST_CASE("std::pair<std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector2d>> adaptedHermiteInterpolant" *
      doctest::description("adaptedHermiteInterpolant()")) {
    // clang-format on
    auto c = [](double t) -> Eigen::Vector2d {
      Eigen::Vector2d ct;
      ct << std::cos(2 * PI * t) + 2. / 3. * std::cos(4 * PI * t),
          3. / 2. * std::sin(2 * PI * t);
      return ct;
    };
    Eigen::VectorXd x = Eigen::VectorXd::Constant(5, 0.5);
    std::pair<std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector2d>> sol =
        adaptedHermiteInterpolant(c, 5, data.x, 1e-3);
    std::pair<std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector2d>> stud =
        adaptedHermiteInterpolant_TEST(c, 5, data.x, 1e-3);
    REQUIRE(sol.first.size() == stud.first.size());
    REQUIRE(sol.second.size() == stud.second.size());
    double err = 0;
    for (unsigned int i = 0; i < sol.first.size(); i++) {
      err += (sol.first[i] - stud.first[i]).norm();
    }
    CHECK(err == doctest::Approx(0.).epsilon(1e-9));
  }

  TEST_CASE("void plotKite" * doctest::description("skipped") *
            doctest::skip()) {}
}
