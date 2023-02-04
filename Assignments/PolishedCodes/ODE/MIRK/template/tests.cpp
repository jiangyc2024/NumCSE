#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <Eigen/Core>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("MIRK") {
  TEST_CASE("Eigen::VectorXd Newton2Steps" *
            doctest::description("Newton method")) {
    Eigen::Vector2d x(-1.0, 1.0);

    // f(x,y) = [e^x,e^y]
    auto f = [](Eigen::Vector2d x) {
      return Eigen::Vector2d(std::exp(x(0)), std::exp(x(1)));
    };

    auto df = [](Eigen::Vector2d x) {
      Eigen::Matrix2d df;
      df << std::exp(x(0)), 0, 0, std::exp(x(1));
      return df;
    };

    // For exp(y): f(y) =  f'(y) for all y and one Newton step satisfies
    // y_{n+1} = y_n - 1
    Eigen::Vector2d x_2ref(-3.0, -1.0);

    // Compute error
    Eigen::Vector2d x_2 = Newton2Steps_TEST(f, df, x);
    REQUIRE(x_2.size() == x_2ref.size());
    CHECK((x_2 - x_2ref).norm() == doctest::Approx(0.).epsilon(1e-7));
  }

  TEST_CASE("double MIRKStep" * doctest::description("step")) {
    constexpr double lam = -50;
    auto f = [](double x) -> double { return lam * x; };
    auto df = [](double x) -> double { return lam; };
    constexpr double y0 = 1.;
    constexpr double h = 0.001;
    auto sol = [](double x) -> double { return y0 * std::exp(lam * x); };

    const double y1 = MIRKStep_TEST(f, df, y0, h);
    CHECK(y1 == doctest::Approx(sol(h)).epsilon(1e-5));
  }

  TEST_CASE("double MIRKSolve" * doctest::description("whole solve")) {
    constexpr double lam = -50;
    auto f = [](double x) -> double { return lam * x; };
    auto df = [](double x) -> double { return lam; };
    constexpr double y0 = 1.;
    auto sol = [](double x) -> double { return y0 * std::exp(lam * x); };
    constexpr double T = 0.1;
    constexpr unsigned int M = 100;

    const double yT = MIRKSolve_TEST(f, df, y0, T, M);
    CHECK(yT == doctest::Approx(sol(T)).epsilon(1e-5));
  }

  TEST_CASE("void cvgMIRK" * doctest::description("convergence")) {
    MESSAGE("This function isn't tested. Run the program to see its output.");
  }
}