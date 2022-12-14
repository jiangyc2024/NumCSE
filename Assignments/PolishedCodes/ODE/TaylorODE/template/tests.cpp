#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <iomanip>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("TaylorIntegrator") {
  TEST_CASE("std::vector<Eigen::Vector2d> SolvePredPreyTaylor" *
            doctest::description("SolvePredPreyTaylor()")) {
    constexpr double alpha1(2), beta1(0.5), alpha2(3), beta2(1), T(10);
    Eigen::Vector2d y0(2, 1);
    constexpr unsigned int M = 100;
    auto sol = SolvePredPreyTaylor(alpha1, beta1, alpha2, beta2, T, y0, M);
    auto stud =
        SolvePredPreyTaylor_TEST(alpha1, beta1, alpha2, beta2, T, y0, M);

    CHECK((sol.back() - stud.back()).norm() ==
          doctest::Approx(0).epsilon(1e-6));
  }

  TEST_CASE("void PrintErrorTable" *
            doctest::description("helper for printing") * doctest::skip()) {}

  TEST_CASE("double testCvgTaylorMethod" *
            doctest::description(
                "TestCvgTaylorMethod() returns the convergence rate")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
