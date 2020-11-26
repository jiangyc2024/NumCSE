#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

auto f = [](double x) { return 1. / (1. + x * x); };
VectorXd t = VectorXd::Random(10);
CubicHermiteInterpolant s_zero(f, t, Slope::Zero);
CubicHermiteInterpolant s_reconstr(f, t, Slope::Reconstructed);
CubicHermiteInterpolant_TEST s_zero_TEST(f, t, Slope::Zero);
CubicHermiteInterpolant_TEST s_reconstr_TEST(f, t, Slope::Reconstructed);

TEST_SUITE("PCHI") {
  TEST_CASE("VectorXd CubicHermiteInterpolant::eval const " * doctest::description("evaluation function")) {
    VectorXd x = VectorXd::Random(10);
    VectorXd sol = s_zero.eval(x);
    VectorXd stud = s_zero_TEST.eval(x);

    CHECK((stud - sol).lpNorm<Infinity>() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::VectorXd fppchip" * doctest::description("...") *
            doctest::skip()) {}

  TEST_CASE("std::vector<double> fppchipConvergence" *
            doctest::description("convergence of flat slope")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("std::vector<double> rspchipConververgence" *
            doctest::description("convergence of reconstructed slope")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
