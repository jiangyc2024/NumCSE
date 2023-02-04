#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

// includes for test data
#include "clenshawcurtisfejer.hpp"
#include <Eigen/Dense>
#include <cmath>

TEST_SUITE("Clenshaw Curtis Fejer Quadrature Rule") {

  TEST_CASE("CCF_QuadRule_Fast()" *
            doctest::description(
                "Checking weights and nodes, not efficiency. This test "
                "requires CCFQuadRule class to be unmodified.")) {
    constexpr unsigned int n = 16;
    const CCFQuadRule_Fast ccfqr_fast(n);
    const CCFQuadRule ccfqr(n);
    CHECK((ccfqr_fast.weights() - ccfqr.weights()).norm() ==
          doctest::Approx(0.).epsilon(1e-10));
    CHECK((ccfqr_fast.nodes() - ccfqr.nodes()).norm() ==
          doctest::Approx(0.).epsilon(1e-10));
  }
}
