#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    n = 10;

    f = [](double t) { return 1 / (1 + std::pow(t, 2)); };
  }

  unsigned int n;
  std::function<double(double)> f;
};

TestData data;

TEST_SUITE("ImproperIntegrals") {
  TEST_CASE("double quadinf" * doctest::description("quadinf")) {
    const double sol = quadinf(data.n, data.f);
    const double stud = quadinf_TEST(data.n, data.f);

    CHECK(sol == doctest::Approx(stud).epsilon(1e-8));
  }

  TEST_CASE("double quad" * doctest::description("Optional helper function") * doctest::skip()) {}

  TEST_CASE("void cvgQuadInf" * doctest::description("Convergence")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
