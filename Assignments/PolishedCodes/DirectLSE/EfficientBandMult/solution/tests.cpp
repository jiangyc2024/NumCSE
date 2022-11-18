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
  TEST_CASE("Eigen::VectorXd multAx" * doctest::description("A*x")) {
    Eigen::VectorXd x = Eigen::VectorXd::Random(data.n);
    Eigen::VectorXd sol, stud;

    sol = multAx(data.a, data.b, x);
    stud = multAx_TEST(data.a, data.b, x);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("Eigen::VectorXd solvelseAupper" *
            doctest::description("A*x = r with b = 0")) {
    Eigen::VectorXd r = Eigen::VectorXd::Random(data.n);
    Eigen::VectorXd sol, stud;

    sol = solvelseAupper(data.a, r);
    stud = solvelseAupper_TEST(data.a, r);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("Eigen::VectorXd solvelseA" * doctest::description("A*x = r")) {
    Eigen::VectorXd r = Eigen::VectorXd::Random(data.n);
    Eigen::VectorXd sol, stud;

    sol = solvelseA(data.a, data.b, r);
    stud = solvelseA_TEST(data.a, data.b, r);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("Eigen::VectorXd solvelseASparse" *
            doctest::description("A*x = r using Eigen sparse methods")) {
    Eigen::VectorXd r = Eigen::VectorXd::Random(data.n);
    Eigen::VectorXd sol, stud;

    sol = solvelseASparse(data.a, data.b, r);
    stud = solvelseASparse_TEST(data.a, data.b, r);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }
}