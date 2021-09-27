#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
    M = 100;
    xa = Eigen::VectorXd::Random(M);
  }

  unsigned int M;
  Eigen::VectorXd xa;
};

TestData data;
constexpr double eps = 1e-10;

TEST_SUITE("StructuredMatrixVector") {
  TEST_CASE("void multAminSlow" * doctest::description("skipped") *
            doctest::skip()) {}

  TEST_CASE("void multAmin" * doctest::description("multAmin")) {
    Eigen::VectorXd sol, stud;

    multAmin(data.xa, sol);
    multAmin_TEST(data.xa, stud);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("void multAmin_runtime" *
            doctest::description("multAmin_runtime")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("Eigen::MatrixXd multABunitv" *
            doctest::description("multABunitv")) {
    Eigen::MatrixXd sol, stud;

    sol = multABunitv();
    stud = multABunitv_TEST();

    REQUIRE(sol.rows() == stud.rows());
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }
}
