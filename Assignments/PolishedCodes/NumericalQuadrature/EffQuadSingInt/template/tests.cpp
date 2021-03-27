#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <cmath>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("EffQuadSingInt") {
  TEST_CASE("double quadsingint" * doctest::description("(8-9.c)")) {
    // Test on monomials
    int n = 0;  // monomial degree
    std::vector<double> ex = {M_PI_2,     0., M_PI_2 / 4,     0,
                              M_PI_2 / 8, 0,  M_PI_2 * 5 / 64};
    auto f = [&n](double t) { return std::pow(t, n); };

    // do not test against solution but against exact values
    for (n = 0; n < ex.size(); ++n) {
      CHECK(quadsingint_TEST(f, 25) == doctest::Approx(ex[n]).epsilon(1e-10));
    }
  }

  TEST_CASE("void tabAndPlotQuadErr" * doctest::description("(8-9.e)")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
