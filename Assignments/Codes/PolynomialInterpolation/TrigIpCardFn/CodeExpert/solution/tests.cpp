#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() { max_pt = 6; }
  unsigned int max_pt;
};

TestData data;

TEST_SUITE("TrigIpCardFn") {
  TEST_CASE("void plot_basis" * doctest::description("plot basis functions")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("void plot_lam" *
            doctest::description("plot lebesgur coefficients")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE(
      "void trigipequid" *
      doctest::description(
          "efficient fft-based computation of coefficients in expansion")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("void trigpolyvalequid" *
            doctest::description("evaluation of trigonometric interpolation "
                                 "polynomial in equidistant points")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("double trigIpL" *
            doctest::description("calculate the lebesgue constants")) {
    for (unsigned int i = 1 << 2; i < (1 << data.max_pt); i = i << 1) {
      double sol = trigIpL(i);
      double stud = trigIpL_TEST(i);

      CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
    }
  }
}
