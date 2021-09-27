#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("NonLinElectr") {
  TEST_CASE("void circuit" * doctest::description("circuit")) {
    constexpr unsigned int n = 20;
    constexpr double alpha = 8;
    constexpr double beta = 1;
    Eigen::VectorXd Uin = Eigen::VectorXd::LinSpaced(n, 0, 20);
    Eigen::VectorXd Uout_stud(n), Uout_reference(n);
    circuit_TEST(alpha, beta, Uin, Uout_stud);
    circuit(alpha, beta, Uin, Uout_reference);
    
    CHECK((Uout_stud - Uout_reference).norm() == doctest::Approx(0.).epsilon(1e-9));
  }
  
  TEST_CASE("void plotU" * doctest::description("plotting")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}

