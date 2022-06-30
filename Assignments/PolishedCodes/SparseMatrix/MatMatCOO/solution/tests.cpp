#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Core>
#include <Eigen/SparseCore>

using Trip = Eigen::Triplet<double>;
using TripVec = std::vector<Trip>;

TEST_SUITE("MatMatCOO") {
  TEST_CASE("TripVec Mat2COO" * doctest::description("conversion")) {}

  TEST_CASE("TripVec COOprod_naive" * doctest::description("naive product")) {}

  TEST_CASE("TripVec COOprod_effic" *
            doctest::description("efficient product")) {}
}
