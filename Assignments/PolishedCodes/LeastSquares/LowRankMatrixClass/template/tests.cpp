#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"
#include <iostream>

// includes for test data
#include <Eigen/Dense>

struct TestData {
  Eigen::Matrix<double, 4, 2> A;
  Eigen::Matrix<double, 3, 2> B;
  Eigen::Matrix<double, 3, 2> C;
  Eigen::Matrix<double, 4, 2> D;
  Eigen::Matrix<double, 3, 1> E;
  Eigen::Matrix<double, 4, 1> F;
  Eigen::Matrix<double, 3, 4> G;
  
  double tol;

  TestData() {
    tol = 1e-4;

    A << 1,2,3,4,5,6,7,8;
    B << 9,0,1,2,3,4;
    C << 5,6,7,8,9,0;
    D << 179,94,437,250,695,406,953,562;
    E << 9.+1e-9,1+1e-9,3+1e-9;
    F << -1,-3,-5,-7;  
    G << 0,0,0,0,4,8,12,16,8,16,24,32;
  }
} data;

struct MasterSolution {
  MatrixLowRank L;
  MatrixLowRank M;
  MatrixLowRank N;

  MasterSolution() {
    L = MatrixLowRank(data.B,data.A);
    M = MatrixLowRank(data.A,data.B);
    N = MatrixLowRank(data.E,data.F);
  }
} sol;

struct StudentSolution {
  MatrixLowRank_TEST L;
  MatrixLowRank_TEST M;
  MatrixLowRank_TEST N;

  StudentSolution() {
    L = MatrixLowRank_TEST(data.B,data.A);
    M = MatrixLowRank_TEST(data.A,data.B);
    N = MatrixLowRank_TEST(data.E,data.F);
  }
} stud;


TEST_SUITE("LowRankMatrix") {
  // Test * operator
  TEST_CASE("Eigen::MatrixXd MatrixLowRank_TEST::operator*()"
      * doctest::description("operator*()")) {
    
    Eigen::MatrixXd MC_sol = sol.M * data.C;
    Eigen::MatrixXd MC_stud = stud.M * data.C;

    double err_sol = (MC_sol - data.D).norm();
    double err_stud = (MC_stud - data.D).norm();

    std::cout << "err_sol=" << err_sol << " err_stud" << err_stud << "\n\n";

    std::cout << "err_sol=" << err_sol << " err_stud=" << err_stud <<"\n\n";

    REQUIRE(MC_stud.cols() == data.D.cols());
    REQUIRE(MC_stud.rows() == data.D.rows());
    //CHECK( (MC - data.D).norm() < doctest::Approx(0.0).epsilon(data.tol) );
    CHECK ( std::abs(err_sol - err_stud) == doctest::Approx(0.0).epsilon(data.tol) );
  }
  
  // Test *= operator
  TEST_CASE("Eigen::MatrixXd MatrixLowRank_TEST::operator*=()"
      * doctest::description("operator*=()")) {
    
    sol.M *= data.C;
    stud.M *= data.C;

    // Multiply with identity so we actualy get an eigen matrix to work with.
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(sol.M.cols(), sol.M.cols());
    Eigen::MatrixXd Me_sol = sol.M * Id; // Eigen version of M.

    // Multiply with identity so we actualy get an eigen matrix to work with.
    Id = Eigen::MatrixXd::Identity(stud.M.cols(), stud.M.cols());
    Eigen::MatrixXd Me_stud = stud.M * Id; // Eigen version of M.
    
    double err_sol = (Me_sol - data.D).norm();
    double err_stud = (Me_stud - data.D).norm();

    REQUIRE(Me_stud.cols() == data.D.cols());
    REQUIRE(Me_stud.rows() == data.D.rows());
    CHECK( std::abs(err_sol - err_stud) == doctest::Approx(0.0).epsilon(data.tol) );
  }

  // Test addTo()
  TEST_CASE("MatrixLowRank_TEST& MatrixLowRank_TEST::addTo"
      * doctest::description("addTo()")) {
 
    sol.N.addTo(sol.L);
    stud.N.addTo(stud.L);
    
    // Multiply with identity so we actualy get an eigen matrix to work with.
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(sol.N.cols(), sol.N.cols());
    Eigen::MatrixXd Ne_sol = sol.N * Id; // Eigen version of M.
    
    // Multiply with identity so we actualy get an eigen matrix to work with.
    Id = Eigen::MatrixXd::Identity(stud.N.cols(), stud.N.cols());
    Eigen::MatrixXd Ne_stud = stud.N * Id; // Eigen version of M.
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(data.G, Eigen::ComputeFullU | Eigen::ComputeFullV);
    svd.setThreshold(1e-6);

    double err_sol = (Ne_sol - data.G).norm();
    double err_stud = (Ne_stud - data.G).norm();


    REQUIRE(Ne_stud.cols() == data.G.cols());
    REQUIRE(Ne_stud.rows() == data.G.rows());
    REQUIRE(stud.N.rank() == svd.rank());
    REQUIRE(stud.N.cols() == data.G.cols());
    REQUIRE(stud.N.rows() == data.G.rows());

    CHECK( std::abs(err_sol - err_stud)  == doctest::Approx(0.0).epsilon(data.tol) );
  }
}
