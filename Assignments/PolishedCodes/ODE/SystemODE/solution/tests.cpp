#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    y0.resize(3);
    y1.resize(3);
    y0 << -1, 1, 2;

    h = 0.5;

    f = [](Eigen::VectorXd y) {
      Eigen::VectorXd fy(3);
      fy << y(1) * y(2), y(0) * y(1), 3 * y(2);
      return fy;
    };
  }
  Eigen::VectorXd y0;
  Eigen::VectorXd y1;

  double h;

  std::function<Eigen::VectorXd(Eigen::VectorXd)> f;
};

TestData data;

TEST_SUITE("SystemODE") {
  TEST_CASE("void rk4step" * doctest::description("Step")) {
    Eigen::VectorXd sol;
    Eigen::VectorXd stud;

    auto f_copy = data.f;

    rk4step(std::move(data.f), data.h, data.y0, sol);
    rk4step_TEST(std::move(f_copy), data.h, data.y0, stud);

    bool samesize = sol.size() == stud.size();
    REQUIRE(samesize);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("double testcvgRK4" * doctest::description("Convergence rate") *
            doctest::skip()) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
