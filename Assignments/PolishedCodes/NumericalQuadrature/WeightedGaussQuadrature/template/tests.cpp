#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <cmath>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("Weighted Gauss quadrature") {
  TEST_CASE("double quadU" * doctest::description("Integration")) {
    auto f = [](double x) { return std::pow(x, 2); };

    // check for exactness where it is known
    CHECK(quadU_TEST(f, 3) == doctest::Approx(M_PI / 8.).epsilon(1e-9));

    auto g = [](double t) { return 1. / (2. + std::exp(3. * t)); };
    CHECK(quadU_TEST(g, 18) == doctest::Approx(quadU(g, 18)).epsilon(1e-9));
  }

  TEST_CASE("void testQuadU" *
            doctest::description("Plotting and tabulating")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
