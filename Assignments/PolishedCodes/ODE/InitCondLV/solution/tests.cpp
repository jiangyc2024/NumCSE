#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    u0 = 2.8;
    v0 = 1.5;
    T = 2;
  }

  double u0;
  double v0;
  double T;
};

TestData data;

TEST_SUITE("InitCondLV") {
  TEST_CASE("std::pair<Eigen::Vector2d, Eigen::Matrix2d> PhiAndW" *
            doctest::description("Test of PhiAndW")) {
    auto [phi_sol, W_sol] = PhiAndW(data.u0, data.v0, data.T);
    auto [phi_stud, W_stud] = PhiAndW_TEST(data.u0, data.v0, data.T);

    CHECK((phi_sol - phi_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    CHECK((W_sol - W_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("std::pair<double, double> findInitCond" *
            doctest::description("findInitCond")) {
    auto [u0_sol, v0_sol] = findInitCond();
    auto [u0_stud, v0_stud] = findInitCond_TEST();

    CHECK(u0_sol == doctest::Approx(u0_stud).epsilon(1e-6));
    CHECK(v0_sol == doctest::Approx(v0_stud).epsilon(1e-6));
  }
}
