#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>
#include <cmath>

#include "arccosquad.hpp"

TEST_SUITE("Smooth integrand by transformation") {
  TEST_CASE("void testConvGaussQuad" *
            doctest::description("Convergence table for gaussquad.")) {
    MESSAGE(
        "This function simply prints a table, and thus is not tested. Run "
        "the program to see its output.");
  }

  TEST_CASE("double arccosWeightedQuad" *
            doctest::description("Approximates I(f) with exponential "
                                 "convergence so a small n should suffice.")) {
    // with f(t) = t^3
    unsigned int n = 16;
    auto f = [](double t) { return t * t * t; };
    double sol = -5 * M_PI / 32.0;
    double ans = arccosWeightedQuad(f, n);

    CHECK(std::abs(sol - ans) == doctest::Approx(0.).epsilon(1e-10));

    // with g(t) = 1/(t - 2)
    n = 30;
    auto g = [](double t) { return 1. / (t - 2.); };
    sol = -1.491633131460735;
    ans = arccosWeightedQuad(g, n);

    CHECK(std::abs(sol - ans) == doctest::Approx(0.).epsilon(1e-10));
  }

  TEST_CASE("void testConvTrfGaussQuad" *
            doctest::description("Convergence table for arccosWeightedQuad.")) {
    MESSAGE(
        "This function simply prints a table, and thus is not tested. Run "
        "the program to see its output.");
  }
}
