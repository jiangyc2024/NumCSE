#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <random>

#include "auxiliary_functions.hpp"
#include "doctest.h"
#include "gaussseidelcrs.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
    // Initialising mat1, Mat1 and b1
    b1.resize(3);
    mat1 = {3, 3, std::vector<double>{1, 2, 3, 4},
            std::vector<unsigned>{0, 2, 1, 2},
            std::vector<unsigned>{0, 2, 3, 4}};
    Mat1.resize(3, 3);
    Mat1 << 1, 0, 2, 0, 3, 0, 0, 0, 4;
    b1 << 1, 1, 1;

    // Initializing mat2, Mat2 and b2
    mat2 = {4, 4, std::vector<double>{1, 2, 3, 4},
            std::vector<unsigned>{0, 1, 2, 3},
            std::vector<unsigned>{0, 1, 2, 3, 4}};
    b2.resize(4);
    Mat2.resize(4, 4);
    Mat2 << 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 4;
    b2 << 1, 1, 1, 1;

    // Initializing mat3
    mat3 = {4, 4, std::vector<double>{1, 2, 3, 4},
            std::vector<unsigned>{1, 1, 2, 3},
            std::vector<unsigned>{0, 1, 2, 3, 4}};
  }

  // Declare variables here
  CRSMatrix mat1;
  Eigen::MatrixXd Mat1;
  Eigen::VectorXd b1;

  CRSMatrix mat2;
  Eigen::MatrixXd Mat2;
  Eigen::VectorXd b2;

  CRSMatrix mat3;
};

TestData data;

TEST_SUITE("Gauss Seidel CRS") {
  TEST_CASE("GaussSeidelstep_crs()" *
            doctest::description("Testing the output")) {
    // Test case I
    Eigen::VectorXd x = Eigen::VectorXd::Random(3);
    Eigen::VectorXd X(x);
    bool pass1 = GaussSeidelstep_crs(data.mat1, data.b1, x);
    CHECK(pass1);

    // Obtaining the result from a dense matrix implementation
    bool pass_dense1 = GaussSeidelstep_generic(data.Mat1, data.b1, X);
    REQUIRE((X - x).norm() == doctest::Approx(0.).epsilon(1e-13));

    // Test case II
    Eigen::VectorXd y = Eigen::VectorXd::Random(4);
    Eigen::VectorXd Y(y);
    bool pass2 = GaussSeidelstep_crs(data.mat2, data.b2, y);
    CHECK(pass2);

    // Obtaining the result from a dense matrix implementation
    bool pass_dense2 = GaussSeidelstep_generic(data.Mat2, data.b2, Y);
    REQUIRE((Y - y).norm() == doctest::Approx(0.).epsilon(1e-13));

    // Test case II
    Eigen::VectorXd z = Eigen::VectorXd::Random(4);
    bool pass3 = GaussSeidelstep_crs(data.mat3, data.b2, z);
    REQUIRE(!pass3);
  }

  TEST_CASE("GaussSeidel_iteration()" * doctest::description("DESCRIPTION")) {
    // Using the student implementation
    Eigen::VectorXd x = Eigen::VectorXd::Random(3);
    bool pass = GaussSeidel_iteration(data.mat1, data.b1, x);
    Eigen::Vector3d sol1(1. / 2, 1. / 3, 1. / 4);

    CHECK(pass);
    REQUIRE((x - sol1).norm() == doctest::Approx(0.).epsilon(1e-6));

    // Using the student implementation
    Eigen::VectorXd y = Eigen::VectorXd::Random(4);
    bool pass_rand = GaussSeidel_iteration(data.mat2, data.b2, y);
    Eigen::VectorXd sol2 = data.Mat2.lu().solve(data.b2);

    CHECK(pass_rand);
    REQUIRE((sol2 - y).norm() == doctest::Approx(0.).epsilon(1e-7));
  }
}
