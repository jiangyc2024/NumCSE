#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
    z = Eigen::VectorXd::Random(15);
    g = Eigen::VectorXd::Random(15);
  }

  Eigen::VectorXd z, g;
};

TestData data;

TEST_SUITE("LeastSquaresFrobenius") {
  TEST_CASE("Eigen::MatrixXd min_frob" *
            doctest::description("augmented least squares")) {
    Eigen::MatrixXd sol = min_frob(data.z, data.g);
    Eigen::MatrixXd stud = min_frob_TEST(data.z, data.g);

    REQUIRE(sol.rows() == stud.rows());
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("bool testMformula" * doctest::description("skipped") *
            doctest::skip()) {}
}
