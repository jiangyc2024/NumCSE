#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    a = Eigen::VectorXd::Random(n - 1);
    b = Eigen::VectorXd::Random(n - 2);
  }
  const unsigned int n = 30;
  Eigen::VectorXd a, b;
};

TestData data;
constexpr double eps = 1e-12;

TEST_SUITE("EfficientBandMult") {
  TEST_CASE("void multAx" * doctest::description("A*x")) {
    Eigen::VectorXd x = Eigen::VectorXd::Random(data.n);
    Eigen::VectorXd sol, stud;

    multAx(data.a, data.b, x, sol);
    multAx_TEST(data.a, data.b, x, stud);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("void solvelseAupper" *
            doctest::description("A*x = r with b = 0")) {
    Eigen::VectorXd r = Eigen::VectorXd::Random(data.n);
    Eigen::VectorXd sol, stud;

    solvelseAupper(data.a, r, sol);
    solvelseAupper_TEST(data.a, r, stud);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("void solvelseA" * doctest::description("A*x = r")) {
    Eigen::VectorXd r = Eigen::VectorXd::Random(data.n);
    Eigen::VectorXd sol, stud;

    solvelseA(data.a, data.b, r, sol);
    solvelseA_TEST(data.a, data.b, r, stud);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("void solvelseAEigen" *
            doctest::description("A*x = r using Eigen")) {
    Eigen::VectorXd r = Eigen::VectorXd::Random(data.n);
    Eigen::VectorXd sol, stud;

    solvelseAEigen(data.a, data.b, r, sol);
    solvelseAEigen_TEST(data.a, data.b, r, stud);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }
}