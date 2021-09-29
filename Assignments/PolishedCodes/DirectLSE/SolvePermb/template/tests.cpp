#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
    A = Eigen::MatrixXd::Random(n, n);
    b = Eigen::VectorXd::Random(n);
    solvpermb_on3(A, b, X_sol);
  }
  Eigen::MatrixXd A, X_sol;
  Eigen::VectorXd b;
  const unsigned int n = 25;
};

TestData data;
constexpr double eps = 1e-9;

TEST_SUITE("SolvePermb") {
  TEST_CASE("void shift" *
            doctest::description("Helper function - not tested") *
            doctest::skip()) {}

  TEST_CASE("void solvpermb" * doctest::description("Naive implementation")) {
    Eigen::VectorXd b = data.b;
    Eigen::MatrixXd X_stud;
    solvpermb_TEST(data.A, b, X_stud);
    REQUIRE(X_stud.rows() == data.X_sol.rows());
    REQUIRE(X_stud.cols() == data.X_sol.cols());
    CHECK((X_stud - data.X_sol).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("void solvpermb_on3" *
            doctest::description("O(n^3) implementation")) {
    Eigen::VectorXd b = data.b;
    Eigen::MatrixXd X_stud;
    solvpermb_on3_TEST(data.A, b, X_stud);
    REQUIRE(X_stud.rows() == data.X_sol.rows());
    REQUIRE(X_stud.cols() == data.X_sol.cols());
    CHECK((X_stud - data.X_sol).norm() == doctest::Approx(0.).epsilon(eps));
  }
}
