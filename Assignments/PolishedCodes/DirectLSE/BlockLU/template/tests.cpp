#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

constexpr double eps = 1e-10;

TEST_SUITE("BlockLU") {
  TEST_CASE("Eigen::VectorXd solve_R" * doctest::description("optional") *
            doctest::skip()) {}

  TEST_CASE("Eigen::VectorXd solve_LSE" * doctest::description("solve_LSE()")) {
    constexpr unsigned int n = 10;
    Eigen::VectorXd v, u, b;
    u = v = Eigen::VectorXd::Random(n);
    b = Eigen::VectorXd::Random(n + 1);
    // upper triangular matrix
    Eigen::MatrixXd R(n, n);
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = i; j < n; ++j) {
        R(i, j) = rand();  // Bad RNG, but sufficient here
      }
    }
    R /= RAND_MAX;  // "norm" R for numerical stability
                    // Build matrix A for Eigensolver
    Eigen::MatrixXd A(n + 1, n + 1);
    A << R, v, u.transpose(), 0;

    const double error =
        (solve_LSE_TEST(R, v, u, b) - A.colPivHouseholderQr().solve(b)).norm();

    CHECK(error == doctest::Approx(0.).epsilon(eps));
  }
}
