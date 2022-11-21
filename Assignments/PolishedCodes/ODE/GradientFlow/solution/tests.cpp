#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

using namespace Eigen;

struct TestData {
  TestData() {
    h = .1;
    n = 30;
    Y0 << 10, 1;
  }

  double h;
  unsigned int n;
  Eigen::Vector2d Y0;
};

TestData data;

// Test function f:
Eigen::Vector2d f(const Vector2d &x) {
  Vector2d y;
  y << x(0)*x(0)+x(1)*x(1), x(0)*x(1);
  return y;
}

// Jacobian of f:
Eigen::Matrix2d J(const Vector2d &x) {
  Matrix2d Df;
  Df << 2*x(0), 2*x(1), x(1), x(0);
  return Df;
}

TEST_SUITE("GradientFlow") {
  TEST_CASE("std::array<Eigen::VectorXd, 5> computeStages" *
            doctest::description("Check computed stages")) {
    auto sols = computeStages(f,J,data.Y0, data.h);
    auto studs = computeStages_TEST(f,J,data.Y0, data.h);

    for (unsigned i = 1 ; i < 5 ; ++i) {
      auto sol = sols[i];
      auto stud = studs[i];
      bool samesize = sol.rows() == stud.rows() && sol.cols() == stud.cols();
      CHECK(samesize);
      if (samesize) {
        CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
      }
    }
  }

  TEST_CASE("Eigen::VectorXd discEvolSDIRK" *
            doctest::description(
                "Checking output of the discrete evolution operator")) {
    auto sol = discEvolSDIRK(f,J,data.Y0, data.h);
    auto stud = discEvolSDIRK_TEST(f,J,data.Y0, data.h);

    bool samesize = sol.rows() == stud.rows() && sol.cols() == stud.cols();
    CHECK(samesize);
    if (samesize) {
      CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("std::vector<Eigen::VectorXd> solveGradientFlow" *
            doctest::description(
                "Checking last entry")) {
                  Eigen::VectorXd d = Eigen::VectorXd::Random(5,1);
                  Eigen::VectorXd y = Eigen::VectorXd::Random(5,1);
                  double lambda = 0.5;
                  double T = 10;
                  unsigned M = 50;
    auto sols = solveGradientFlow(d,lambda,y,T,M);
    auto studs = solveGradientFlow_TEST(d,lambda,y,T,M);
    auto sol = sols.back();
    auto stud = studs.back();

    bool samesize = sol.rows() == stud.rows() && sol.cols() == stud.cols();
    CHECK(samesize);
    if (samesize) {
      CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("Eigen::VectorXd solveGenStageEquation" *
            doctest::description("solveGenStageEquation()")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("Eigen::MatrixXd ButcherMatrix" *
            doctest::description("ButcherMatrix")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
