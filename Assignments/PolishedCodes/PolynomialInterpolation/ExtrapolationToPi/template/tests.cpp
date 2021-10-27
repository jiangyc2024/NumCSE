#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

struct TestData2 {
  Eigen::Matrix<unsigned int, 6, 1> K;

  TestData2() {
    K << 2, 3, 5, 10, 15, 20;
  }
} data;

TEST_SUITE("ExtrapolationToPi") {
  TEST_CASE("double extrapolate_to_pi" *
            doctest::description("Extrapolate to k")) {
    for (unsigned int k : data.K) {
      double sol = extrapolate_to_pi(k);
      double stud = extrapolate_to_pi_TEST(k);

      CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("void plotExtrapolationError" * doctest::description("Plot error") *
            doctest::skip()) {}
}
