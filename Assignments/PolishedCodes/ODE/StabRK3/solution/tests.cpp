#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    T = 1;
    n = 100;
    Y0 << 15, 1;
  }

  double T;
  unsigned int n;
  Eigen::Vector2d Y0;
};

TestData data;

TEST_SUITE("STABRK3") {
  TEST_CASE("Eigen::Vector2d PredPrey" * doctest::description("Check output")) {
    Eigen::Vector2d sol = PredPrey(data.Y0, data.T, data.n);
    Eigen::Vector2d stud = PredPrey_TEST(data.Y0, data.T, data.n);

    const bool samesize =
        sol.rows() == stud.rows() && sol.cols() == stud.cols();
    REQUIRE(samesize);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("void SimulatePredPrey" *
            doctest::description("SimulatePredPrey()")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
