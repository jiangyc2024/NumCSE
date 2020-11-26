#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"
// #include "pchi.hpp"

auto f = [](double x) { return 1. / (1. + x * x); };
VectorXd t = VectorXd::Random(10);
CubicHermiteInterpolant s_zero(f, t, Slope::Zero);
CubicHermiteInterpolant s_reconstr(f, t, Slope::Reconstructed);
CubicHermiteInterpolant_TEST s_zero_TEST(f, t, Slope::Zero);
CubicHermiteInterpolant_TEST s_reconstr_TEST(f, t, Slope::Reconstructed);

TEST_SUITE("PCHI") {
  TEST_CASE("eval" * doctest::description("evaluation function")) {
    VectorXd x = VectorXd::Random(10);
    VectorXd sol_zero = s_zero.eval(x);
    VectorXd stud_zero = s_zero_TEST.eval(x);
    VectorXd sol_reconstr = s_reconstr.eval(x);
    VectorXd stud_reconstr = s_reconstr_TEST.eval(x);

    CHECK((sol_zero - stud_zero).lpNorm<Infinity>() ==
          doctest::Approx(0.).epsilon(1e-6));
    CHECK((sol_reconstr - stud_reconstr).lpNorm<Infinity>() ==
          doctest::Approx(0.).epsilon(1e-6));
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
