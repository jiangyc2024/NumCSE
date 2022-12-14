#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    T = 1.1;
    p0 << 0.1, 0.0, 0.5;
    q0 << 0.0, 0.1, 0.5;
    M = 5;
    tau = 0.1;
  }
  double T, tau;
  Eigen::Vector3d p0, q0;
  int M;
};

TestData data;

TEST_SUITE("SymplecticTimestepping") {
  TEST_CASE("void sympTimestep" *
            doctest::description("Check the timestepping scheme")) {
    Eigen::Vector2d sol(1, 1);
    Eigen::Vector2d stud(sol);

    sympTimestep(data.tau, sol);
    sympTimestep_TEST(data.tau, stud);

    const bool samesize =
        sol.rows() == stud.rows() && sol.cols() == stud.cols();
    REQUIRE(samesize);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::MatrixXd simulateHamiltonianDynamics" *
            doctest::description("Hamiltonian dynamics timestepper")) {
    auto sol = simulateHamiltonianDynamics(data.p0, data.q0, data.T, data.M);
    auto stud =
        simulateHamiltonianDynamics_TEST(data.p0, data.q0, data.T, data.M);

    const bool samesize =
        sol.rows() == stud.rows() && sol.cols() == stud.cols();
    REQUIRE(samesize);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("void sympTimesteppingODETest" *
            doctest::description("sympTimesteppingODETest()")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
