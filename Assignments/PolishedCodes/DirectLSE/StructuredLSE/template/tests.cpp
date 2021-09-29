#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
    a = Eigen::VectorXd::Random(n);
    b = Eigen::VectorXd::Random(n);
  }
  const unsigned int n = 25;
  Eigen::VectorXd a, b;
};

TestData data;
constexpr double eps = 1e-12;

TEST_SUITE("StructuredLSE") {
  TEST_CASE("Eigen::MatrixXd buildA" *
            doctest::description("Build system matrix")) {
    Eigen::MatrixXd A_sol = buildA(data.a);
    Eigen::MatrixXd A_stud = buildA_TEST(data.a);

    REQUIRE(A_sol.rows() == A_stud.rows());
    REQUIRE(A_sol.cols() == A_stud.cols());
    CHECK((A_sol - A_stud).norm() == doctest::Approx(0.).epsilon(1e-9));
  }

  TEST_CASE("void solveA" * doctest::description("Naive solver")) {
    Eigen::VectorXd x_sol, x_stud;

    solveA(data.a, data.b, x_sol);
    solveA_TEST(data.a, data.b, x_stud);
    REQUIRE(x_sol.size() == x_stud.size());
    CHECK((x_sol - x_stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("void solveA_fast" * doctest::description("Fast solver")) {
    Eigen::VectorXd x_sol, x_stud;

    solveA_fast(data.a, data.b, x_sol);
    solveA_fast_TEST(data.a, data.b, x_stud);
    REQUIRE(x_sol.size() == x_stud.size());
    CHECK((x_sol - x_stud).norm() == doctest::Approx(0.).epsilon(eps));
  }
}
