#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "doctest.h"
#include "xmatmult.hpp"

TEST_SUITE("XMatMult") {
  TEST_CASE("void xmatmult" * doctest::description("xmatmult")) {
    // testing for even n
    unsigned int n = 10;
    Eigen::VectorXd a, y, x(n);
    a = y = Eigen::VectorXd::Random(n, 1);
    // building A for normal Matrix-Vector multiplication O(n*n)
    Eigen::MatrixXd A = a.asDiagonal();
    for (unsigned int i = 0; i < n; ++i) {
      A(n - i - 1, i) = A(i, i);
    }
    xmatmult(a, y, x);

    REQUIRE(x.size() == 10);
    CHECK((x - A * y).norm() == doctest::Approx(0.).epsilon(1e-10));

    // testing for odd n
    n = 11;
    a = y = x = Eigen::VectorXd::Random(n, 1);
    A = a.asDiagonal();
    for (unsigned int i = 0; i < n; ++i) {
      A(n - i - 1, i) = A(i, i);
    }
    xmatmult(a, y, x);

    REQUIRE(x.size() == 11);
    CHECK((x - A * y).norm() == doctest::Approx(0.).epsilon(1e-10));
  }

  TEST_CASE("void compare_times" * doctest::description("timings")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
