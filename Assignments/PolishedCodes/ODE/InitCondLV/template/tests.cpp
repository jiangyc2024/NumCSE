#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

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
    std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW_sol =
        PhiAndW(data.u0, data.v0, data.T);
    std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW_stud =
        PhiAndW_TEST(data.u0, data.v0, data.T);

    Eigen::Vector2d phi_sol = PaW_sol.first;
    Eigen::Vector2d phi_stud = PaW_stud.first;

    CHECK((phi_sol - phi_stud).norm() == doctest::Approx(0.).epsilon(1e-6));

    Eigen::Matrix2d W_sol = PaW_sol.second;
    Eigen::Matrix2d W_stud = PaW_stud.second;

    CHECK((W_sol - W_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::Vector2d findInitCond" *
            doctest::description("findInitCond: y")) {
    Eigen::Vector2d sol = findInitCond();
    Eigen::Vector2d stud = findInitCond_TEST();

    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }
}
