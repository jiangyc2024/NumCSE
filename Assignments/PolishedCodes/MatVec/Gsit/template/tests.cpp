#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "doctest.h"
#include "gsit.hpp"

TEST_SUITE("GSIt") {
  TEST_CASE("void GSIt" * doctest::description("Gauss-Seidel iteration")) {
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
    A = 0.5 * (A + A.transpose()) + 10. * Eigen::MatrixXd::Identity(10, 10);
    const Eigen::VectorXd b = Eigen::VectorXd::Random(10);
    Eigen::VectorXd x = b;
    constexpr double rtol = 1e-8;
    GSIt(A, b, x, rtol);
    REQUIRE(x.size() == 10);
    CHECK((A * x - b).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("double testGSIt" * doctest::description("Test")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
