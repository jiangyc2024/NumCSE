#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <iomanip>
#include <random>

#include "convolutionquadrature.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    // Initialise variables here
    wts = Eigen::VectorXcd::Constant(
        11, std::complex<double>(0.10004885197850513, 0));
  }
  // Declare variables here
  Eigen::VectorXcd wts;
};

TestData data;

TEST_SUITE("Convolution Quadrature") {
  TEST_CASE("Eigen::VectorXcd compute_cq_weights" *
            doctest::description("Testing output. Not testing efficiency")) {
    const double n = static_cast<double>(rand()) / RAND_MAX;
    auto F1 = [&](std::complex<double> x) {
      return std::complex<double>(n, 0);
    };

    auto F2 = [](std::complex<double> x) {
      return std::complex<double>(1, 0) / x;
    };
    constexpr unsigned int N = 10;
    constexpr double tau = .1;

    Eigen::VectorXcd wts1 = compute_cq_weights(F1, N, tau);
    Eigen::VectorXcd ex_wts1 =
        Eigen::VectorXcd::Constant(N + 1, std::complex<double>(0, 0));
    ex_wts1(0) = std::complex<double>(n, 0);

    REQUIRE(wts1.size() == N + 1);
    CHECK((wts1 - ex_wts1).norm() == doctest::Approx(0.).epsilon(1e-7));

    Eigen::VectorXcd wts2 = compute_cq_weights(F2, N, tau, 0.5);
    REQUIRE(wts2.size() == N + 1);
    CHECK((wts2 - data.wts).norm() == doctest::Approx(0.).epsilon(1e-7));
  }

  TEST_CASE("Eigen::VectorXd compute_cq_weights_sqrt" *
            doctest::description("helper") * doctest::skip()) {}
}
