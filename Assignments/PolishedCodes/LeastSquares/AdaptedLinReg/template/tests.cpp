#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    n = 50;
    t = Eigen::VectorXd::LinSpaced(n, 0., 1.);
    noise = Eigen::VectorXd::Random(n);
  }

  std::size_t n;
  Eigen::VectorXd t, noise;
};

TestData data;

TEST_SUITE("AdaptedLinReg") {
  TEST_CASE("Eigen::VectorXd linReg" * doctest::description("linReg()")) {
    Eigen::VectorXd y =
        12 * data.t + Eigen::VectorXd::Constant(data.n, -154) + data.noise;
    Eigen::VectorXd sol, stud;

    sol = linReg(data.t, y);
    stud = linReg_TEST(data.t, y);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("Eigen::VectorXd expFit" * doctest::description("expFit()")) {
    Eigen::VectorXd y =
        17 * Eigen::exp(-3 * data.t.array()).matrix() + 0.5 * data.noise;
    Eigen::VectorXd sol, stud;

    sol = expFit(data.t, y);
    stud = expFit_TEST(data.t, y);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }
}
