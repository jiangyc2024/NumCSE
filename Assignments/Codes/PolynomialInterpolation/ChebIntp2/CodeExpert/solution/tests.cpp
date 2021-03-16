#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    a = 0;
    b = 1;
    N = 100;
    tol = 1e-4;

    f = [](double t) { return std::sin(std::exp(2 * t)); };
  }

  double a;
  double b;
  unsigned int N;
  double tol;

  std::function<double(double)> f;
};

TestData data;

TEST_SUITE("AdaptivePolyIntp") {
  TEST_CASE("VectorXd adaptivepolyintp" * doctest::description("Nodes")) {
    Eigen::VectorXd sol_vec =
        adaptivepolyintp(data.f, data.a, data.b, data.tol, data.N);
    Eigen::VectorXd stud_vec =
        adaptivepolyintp_TEST(data.f, data.a, data.b, data.tol, data.N);

    const bool samesize = (sol_vec.size() == stud_vec.size());
    CHECK(samesize);

    if (samesize) {
      for (std::size_t i = 0; i < sol_vec.size(); i++) {
        const double sol = sol_vec[i];
        const double stud = stud_vec[i];

        CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-9));
      }
    }
  }

  TEST_CASE("void plotInterpolationError" *
            doctest::description("Plot error")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
