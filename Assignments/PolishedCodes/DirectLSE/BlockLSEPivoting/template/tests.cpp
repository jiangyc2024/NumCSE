#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <Eigen/Dense>

#include "blocklsepiv.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    d1 = Eigen::VectorXd::LinSpaced(n, 1, n);
    d2 = -d1;
    c = Eigen::VectorXd::Ones(n);
    b.resize(2 * n);
    b << d1, d1;
    A.resize(2 * n, 2 * n);
    A << static_cast<Eigen::MatrixXd>(d1.asDiagonal()),
        static_cast<Eigen::MatrixXd>(c.asDiagonal()),
        static_cast<Eigen::MatrixXd>(c.asDiagonal()),
        static_cast<Eigen::MatrixXd>(d2.asDiagonal());
  }

  const unsigned int n = 50;
  Eigen::VectorXd d1, d2, c, b;
  Eigen::MatrixXd A;
};

TestData t;
constexpr double eps = 1e-12;

TEST_SUITE("BlockLSEPiv") {
  TEST_CASE("Eigen::VectorXd multA" * doctest::description("multA")) {
    Eigen::VectorXd ye = t.A * t.b;
    Eigen::VectorXd yo = multA(t.d1, t.d2, t.c, t.b);
    CHECK((ye - yo).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("Eigen::VectorXd solveA" * doctest::description("solveA")) {
    Eigen::VectorXd ye = t.A.partialPivLu().solve(t.b);
    Eigen::VectorXd yo = solveA(t.d1, t.d2, t.c, t.b);
    CHECK((ye - yo).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("void numericalExperiment" *
            doctest::description("numerical experiment")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
