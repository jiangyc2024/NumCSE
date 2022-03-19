#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <iomanip>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  Eigen::VectorXd coeff, t, x_vec, y;
  double x;

  TestData() {
    coeff = Eigen::VectorXd::Random(10);
    x = -2.2;
    t = Eigen::VectorXd::LinSpaced(10, -5, 5);
    y = Eigen::VectorXd::Random(10);
    x_vec = Eigen::VectorXd::LinSpaced(20, -10, 10);
  }
} data;

TEST_SUITE("EvaluatingDerivatives") {
  TEST_CASE("std::pair<double, double> evaldp" *
            doctest::description("evaldp()")) {
    auto sol = evaldp(data.coeff, data.x);
    auto stud = evaldp_TEST(data.coeff, data.x);

    CHECK(sol.first == doctest::Approx(stud.first).epsilon(1e-6));
    CHECK(sol.second == doctest::Approx(stud.second).epsilon(1e-6));
  }

  TEST_CASE("std::pair<double, double> evaldp_naive" *
            doctest::description("evaldp_naive()")) {
    auto sol = evaldp_naive(data.coeff, data.x);
    auto stud = evaldp_naive_TEST(data.coeff, data.x);

    CHECK(sol.first == doctest::Approx(stud.first).epsilon(1e-6));
    CHECK(sol.second == doctest::Approx(stud.second).epsilon(1e-6));
  }

  TEST_CASE("bool polyTestTime" *
            doctest::description("polyTestTime() returns true")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("Eigen::VectorXd dipoleval" * doctest::description("dipoleval()")) {
    auto sol = dipoleval(data.t, data.y, data.x_vec);
    auto stud = dipoleval_TEST(data.t, data.y, data.x_vec);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::VectorXd dipoleval_alt" *
            doctest::description("dipoleval_alt()")) {
    auto sol = dipoleval_alt(data.t, data.y, data.x_vec);
    auto stud = dipoleval_alt_TEST(data.t, data.y, data.x_vec);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("bool testDipolEval" *
            doctest::description("testDipolEval() returns true") *
            doctest::skip()) {}

  TEST_CASE("void plotPolyInterpolant" * doctest::description("skipped") *
            doctest::skip()) {}
}
