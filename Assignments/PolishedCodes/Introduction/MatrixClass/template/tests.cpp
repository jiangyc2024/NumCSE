#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

constexpr double eps = 1e-10;

TEST_SUITE("MatrixClass") {
  TEST_CASE("Eigen::Matrix<double, 2, 2> smallTriangular" *
            doctest::description("smallTriangular")) {
    const auto sol = smallTriangular(1, 2, 3);
    const auto stud = smallTriangular_TEST(1, 2, 3);

    REQUIRE(sol.rows() == stud.rows());
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("Eigen::MatrixXd constantTriangular" *
            doctest::description("constantTriangular")) {
    const auto sol = constantTriangular(3, 20);
    const auto stud = constantTriangular_TEST(3, 20);

    REQUIRE(sol.rows() == stud.rows());
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("double casting" * doctest::description("casting")) {
    const auto sol = casting();
    const auto stud = casting_TEST();

    CHECK(sol == doctest::Approx(stud).epsilon(eps));
  }

  TEST_CASE("Eigen::VectorXcd arithmetics" *
            doctest::description("arithmetics")) {
    const auto sol = arithmetics(4);
    const auto stud = arithmetics_TEST(4);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }
}
