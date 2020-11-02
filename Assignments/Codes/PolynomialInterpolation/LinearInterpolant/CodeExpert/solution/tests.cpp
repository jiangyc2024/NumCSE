#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Core>

TEST_SUITE("LinearInterpolant") {
  TEST_CASE("LinearInterpolant [OUT OF CLASS]" *
            doctest::description("Constructor") * doctest::skip()) {}

  TEST_CASE("double operator() [OUT OF CLASS]" *
            doctest::description("evaluation operator")) {
    Eigen::VectorXd t(4);
    t << 2., -1., 4., 1.;
    Eigen::VectorXd y(4);
    y << 3., -1., 4., 2.;

    LinearInterpolant sol(t, y);
    LinearInterpolant_TEST stud(t, y);

    for (double x = -1.; x <= 4; x += 0.05) {
      REQUIRE(sol(x) == doctest::Approx(stud(x)).epsilon(1e-9));
    }
  }
}
