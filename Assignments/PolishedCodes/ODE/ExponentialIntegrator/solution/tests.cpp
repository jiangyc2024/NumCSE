#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    h = .1;
    Y0.resize(1, 1);
    Y0 << 0.7;
  }

  double h;
  Eigen::MatrixXd Y0;
};

TestData data;

double f(const Eigen::VectorXd &y) { return y(0) * (1.0 - y(0)); }
Eigen::MatrixXd df(const Eigen::VectorXd &y) {
  Eigen::MatrixXd dfy(1, 1);
  dfy << 1.0 - 2.0 * y(0);
  return dfy;
}

TEST_SUITE("ExponentialIntegrator") {
  TEST_CASE("Eigen::VectorXd exponentialEulerStep" *
            doctest::description("Check value after step")) {
    Eigen::MatrixXd sol = exponentialEulerStep(data.Y0, f, df, data.h);
    Eigen::MatrixXd stud = exponentialEulerStep_TEST(data.Y0, f, df, data.h);

    const bool samesize =
        sol.rows() == stud.rows() && sol.cols() == stud.cols();
    REQUIRE(samesize);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::MatrixXd phim" * doctest::description("phim()")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("void testExpEulerLogODE" *
            doctest::description("testExpEulerLogODE()")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
