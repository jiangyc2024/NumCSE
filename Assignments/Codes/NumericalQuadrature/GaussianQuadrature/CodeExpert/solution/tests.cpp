#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    N = 50;
    I = 0.870267525725852642;

    f = [](double t) { return std::sinh(t); };
  }

  unsigned int N;
  double I;
  std::function<double(double)> f;
};

TestData data;

TEST_SUITE("GaussionQuadrature") {
  TEST_CASE("double gaussConv" * doctest::description("gaussian quadrature")) {
    const double sol = gaussConv(data.f, data.I, data.N);
    const double stud = gaussConv_TEST(data.f, data.I, data.N);

    CHECK(sol == doctest::Approx(stud).epsilon(1e-8));
  }

  TEST_CASE("double gaussConvCV" *
            doctest::description("gaussian quadrature substitution")) {
    const double sol = gaussConvCV(data.f, data.I, data.N);
    const double stud = gaussConvCV_TEST(data.f, data.I, data.N);

    CHECK(sol == doctest::Approx(stud).epsilon(1e-8));
  }

  TEST_CASE("double integrate" *
            doctest::description("Optional helper function") *
            doctest::skip()) {}
}
