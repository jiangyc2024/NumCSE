#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

TEST_SUITE("NLMatODE") {
  TEST_CASE("double testCvgRKSSM" *
            doctest::description(
                "Compute convergence rate of an RK Single-Step method")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("void cmpCvgRKSSM" *
            doctest::description(
                "Compare convergence rates for RK Single-Step methods")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
