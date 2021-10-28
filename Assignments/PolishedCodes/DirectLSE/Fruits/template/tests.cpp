#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "copy.hpp"
#include "doctest.h"

constexpr double eps = 1e-10;

TEST_SUITE("Fruits") {
  TEST_CASE("Eigen::VectorXd fruitPrice" * doctest::description("fruitPrice")) {
    const auto sol = fruitPrice();
    const auto stud = fruitPrice_TEST();

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }
}
