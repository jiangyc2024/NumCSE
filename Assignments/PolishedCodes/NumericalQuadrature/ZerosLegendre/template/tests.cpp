#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include <iostream>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("ZerosLegendre") {
  TEST_CASE("void legvals" * doctest::description("legvals (8-10.c)")) {
    constexpr unsigned int N = 20;
    constexpr unsigned int n = 10;
    const Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, 0., 10.);
    Eigen::MatrixXd Lx_sol(N, n), Lx_stud(N, n), DLx_sol(N, n), DLx_stud(N, n);
    
    legvals(x, Lx_sol, DLx_sol);
    legvals_TEST(x, Lx_stud, DLx_stud);
    
    CHECK((Lx_sol - Lx_stud).norm() == doctest::Approx(0.).epsilon(1e-9));
    CHECK((DLx_sol - DLx_stud).norm() == doctest::Approx(0.).epsilon(1e-9));
  }
  
  TEST_CASE("double Pnx" * doctest::description("optional helper function") * doctest::skip()) {}
  
  TEST_CASE("Eigen::MatrixXd gaussPts" * doctest::description("gaussPts (8-10.d)")) {
    
    // This case will pass even though the function gives wrong output (in this case for $P_8$). It just checks if the function is implemented like the solution.
    
    constexpr unsigned int n = 8;
    Eigen::MatrixXd sol = gaussPts(n, 1e-10, 1e-12);
    Eigen::MatrixXd stud = gaussPts_TEST(n, 1e-10, 1e-12);
    
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));
  }
  
  TEST_CASE("Eigen::MatrixXd gaussPts_regulaFalsi" * doctest::description("gaussPts with regula falsi (8-10.f)")) {
    constexpr unsigned int n = 8;
    Eigen::MatrixXd sol = gaussPts_regulaFalsi(n, 1e-10, 1e-12);
    Eigen::MatrixXd stud = gaussPts_regulaFalsi_TEST(n, 1e-10, 1e-12);
    
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));
  }
}
