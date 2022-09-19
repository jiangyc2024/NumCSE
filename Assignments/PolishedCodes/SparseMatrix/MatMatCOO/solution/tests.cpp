#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Core>
#include <Eigen/SparseCore>

using Trip = Eigen::Triplet<double>;
using TripVec = std::vector<Trip>;

TEST_SUITE("MatMatCOO") {
  TEST_CASE("Eigen::MatrixXd COO2Mat" *
            doctest::description("conversion to dense") * doctest::skip()) {}

  TEST_CASE("Eigen::MatrixXd randMat" *
            doctest::description("create random binary matrix") *
            doctest::skip()) {}

  TEST_CASE("TripVec Mat2COO" *
            doctest::description("conversion to COO format")) {
    Eigen::MatrixXd randmat = randMat(50, 50, 0.1);
    TripVec randmat_coo = Mat2COO_TEST(randmat);      // convert to COO
    Eigen::MatrixXd randmat_ = COO2Mat(randmat_coo);  // convert back
    REQUIRE(randmat.rows() == randmat_.rows());
    REQUIRE(randmat.cols() == randmat_.cols());
    CHECK((randmat - randmat_).norm() == doctest::Approx(0.).epsilon(1e-7));
  }

  TEST_CASE("TripVec COOprod_naive" * doctest::description("naive product")) {
    Eigen::MatrixXd randA = randMat(100, 50, 0.1);
    Eigen::MatrixXd randB = randMat(50, 10, 0.1);
    Eigen::MatrixXd sol =
        randA * randB;  // normal dense Eigen matrix-matrix product
    TripVec Acoo = Mat2COO(randA);
    TripVec Bcoo = Mat2COO(randB);
    TripVec Ccoo = COOprod_naive_TEST(Acoo, Bcoo);  // product to test
    Eigen::MatrixXd C = COO2Mat(Ccoo);
    REQUIRE(C.rows() == sol.rows());
    REQUIRE(C.cols() == sol.cols());
    CHECK((C - sol).norm() == doctest::Approx(0.).epsilon(1e-7));
  }

  TEST_CASE("TripVec COOprod_effic" *
            doctest::description("efficient product")) {
    Eigen::MatrixXd randA = randMat(100, 50, 0.1);
    Eigen::MatrixXd randB = randMat(50, 10, 0.1);
    Eigen::MatrixXd sol =
        randA * randB;  // normal dense Eigen matrix-matrix product
    TripVec Acoo = Mat2COO(randA);
    TripVec Bcoo = Mat2COO(randB);
    TripVec Ccoo = COOprod_effic_TEST(Acoo, Bcoo);  // product to test
    Eigen::MatrixXd C = COO2Mat(Ccoo);
    REQUIRE(C.rows() == sol.rows());
    REQUIRE(C.cols() == sol.cols());
    CHECK((C - sol).norm() == doctest::Approx(0.).epsilon(1e-7));
  }
}
