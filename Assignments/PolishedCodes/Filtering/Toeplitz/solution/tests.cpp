#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

TEST_SUITE("Toeplitz") {
  TEST_CASE("Eigen::MatrixXd toeplitz" *
            doctest::description("Toeplitz matrix")) {
    Eigen::VectorXd c = Eigen::VectorXd::Random(50);
    Eigen::VectorXd r = Eigen::VectorXd::Random(50);
    c(0) = r(0);

    auto sol = toeplitz(c, r);
    auto stud = toeplitz_TEST(c, r);

    REQUIRE(sol.rows() == stud.rows());
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::VectorXd toepmatmult" * doctest::description("skipped") *
            doctest::skip()) {}

  TEST_CASE("Eigen::VectorXd toepmult" * doctest::description("skipped") *
            doctest::skip()) {}

  TEST_CASE("Eigen::VectorXd ttmatsolve" * doctest::description("skipped") *
            doctest::skip()) {}

  TEST_CASE("Eigen::VectorXd ttrecsolve" * doctest::description("skipped") *
            doctest::skip()) {}

  TEST_CASE("Eigen::VectorXd ttsolve" *
            doctest::description("ttsolve wrapper")) {
    Eigen::VectorXd h = Eigen::VectorXd::LinSpaced(50, 1, 50);
    Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(50, 51, 100);

    auto sol = ttsolve(h, y);
    auto stud = ttsolve_TEST(h, y);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("void runtime_toeplitz" *
            doctest::description("Test and compare runtime")) {
    MESSAGE(
        "This function wasn't tested. Run the program and check the plots.");
  }
  TEST_CASE("using duration_t" * doctest::description("skipped") *
            doctest::skip()) {}
  TEST_CASE("duration_t timing" * doctest::description("skipped") *
            doctest::skip()) {}
  TEST_CASE("void runtime_toeplitz_with_chrono" *
            doctest::description("implemented only in solution") *
            doctest::skip()) {}
}
