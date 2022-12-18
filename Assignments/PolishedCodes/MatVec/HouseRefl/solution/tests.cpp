#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("HouseRefl") {
  TEST_CASE("void houserefl" * doctest::description("Householder reflection")) {
    Eigen::VectorXd v = Eigen::VectorXd::Random(20);
    Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(20, 19);
    houserefl_TEST(v, Z);

    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(19, 19);
    const bool checksize = Z.rows() == 20 && Z.cols() == 19;
    REQUIRE(checksize);
    CHECK((Z.transpose() * Z - Id).norm() == doctest::Approx(0.).epsilon(1e-9));
  }
}