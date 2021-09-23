#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <cmath>
#include <Eigen/Dense>
#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("GaussPts") {
  TEST_CASE("void gaussrule" * doctest::description("not tested; already implemented") * doctest::skip()) {}
  
  TEST_CASE("Eigen::VectorXd comp_g_gausspts" * doctest::description("8-8.b")) {
    auto f = [](double y){ return 1.; };
    const Eigen::VectorXd sol = comp_g_gausspts(f, 20);
    const Eigen::VectorXd stud = comp_g_gausspts_TEST(f, 20);
    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));
  }
  
  TEST_CASE("double testCompGGaussPts" * doctest::description("8-8.c")) {
    CHECK(testCompGGaussPts_TEST() == doctest::Approx(1.).epsilon(1e-12));
  }
}
