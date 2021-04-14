#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("NewtonArctan") {
  TEST_CASE("double newton_arctan" * doctest::description("arctan")) {
    const std::vector<double> X = {1., 1.15, 1.3, 1.6, 2.};
    for (int i = 0; i < X.size(); i++) {
      const double sol = newton_arctan(X[i]);
      const double stud = newton_arctan_TEST(X[i]);

      CHECK(sol == doctest::Approx(stud).epsilon(1e-12));
    }
  }
}
