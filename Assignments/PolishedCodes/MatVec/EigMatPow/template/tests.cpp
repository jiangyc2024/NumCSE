#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "doctest.h"
#include "getit.hpp"

TEST_SUITE("EigMatPow") {
  TEST_CASE("Eigen::VectorXd getit" * doctest::description("y = A^k x")) {
    // Some arbitrary data to test getit
    Eigen::MatrixXd A(4, 4);
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
    Eigen::VectorXd x(4);
    x << 4, 5, 6, 7;
    constexpr unsigned int k = 9;

    // Testing the implementation with some matrix
    Eigen::VectorXd yg = getit(A, x, k);

    // Checking that getit works
    Eigen::VectorXd yp = A.pow(k) * x;
    const double err = (yg - yp).norm() / yp.norm();
    CHECK(err == doctest::Approx(0.).epsilon(1e-7));
  }
}