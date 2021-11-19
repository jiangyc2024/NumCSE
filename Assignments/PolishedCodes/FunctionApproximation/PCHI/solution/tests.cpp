#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

auto f = [](double x) { return 1. / (1. + x * x); };

Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(10,-1,6);
CubicHermiteInterpolant s_reconstr(f, t);
CubicHermiteInterpolant_TEST s_reconstr_TEST(f, t);

TEST_SUITE("PCHI") {
  TEST_CASE("eval" * doctest::description("evaluation function")) {
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(23, -1, 6);
    
    Eigen::VectorXd sol_reconstr = s_reconstr.eval(x);
    Eigen::VectorXd stud_reconstr = s_reconstr_TEST.eval(x);

    CHECK((sol_reconstr - stud_reconstr).lpNorm<Eigen::Infinity>() ==
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
