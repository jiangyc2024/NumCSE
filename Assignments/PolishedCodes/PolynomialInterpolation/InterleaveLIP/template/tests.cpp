#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Core>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    x = Eigen::VectorXd::LinSpaced(11, 0., 10.);
    t.resize(11);
    t[0] = 0;
    t.tail(10).setLinSpaced(10, 0.5, 9.5);
    y.resize(11);
    auto f = [](double x) { return x * x; };
    y = t.unaryExpr(f);
  }

  Eigen::VectorXd x;
  Eigen::VectorXd y;
  Eigen::VectorXd t;
};

TestData data;

TEST_SUITE("InterleaveLIP") {
  TEST_CASE("std::vector<std::size_t> order" *
            doctest::description("Helper function") * doctest::skip()) {}

  TEST_CASE("Eigen::VectorXd tentBasCoeff" *
            doctest::description("tentBasCoeff")) {
    Eigen::VectorXd sol = tentBasCoeff(data.x, data.t, data.y);
    Eigen::VectorXd stud = tentBasCoeff_TEST(data.x, data.t, data.y);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("PwLinIP [OUT OF CLASS]" * doctest::description("Constructor") *
            doctest::skip()) {}

  TEST_CASE("double operator() [OUT OF CLASS]" *
            doctest::description("Evaluation")) {
    PwLinIP sol(data.x, data.t, data.y);
    PwLinIP_TEST stud(data.x, data.t, data.y);
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(1000, 0., 10.);
    Eigen::VectorXd sol_(1000);
    Eigen::VectorXd stud_(1000);
    for (unsigned int i = 0; i < 1000; ++i) {
      sol_[i] = sol(x[i]);
      stud_[i] = stud(x[i]);
    }
    CHECK((sol_ - stud_).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("void plotCardinalBasisFunctions" *
            doctest::description("Plotting")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
