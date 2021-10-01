// NO_CLASS_COPY
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    T.n_rows = m;
    T.n_cols = n;
    T.triplets = {{0, 2, 1}, {0, 4, 1}, {2, 1, 1}, {1, 1, 1},
                  {2, 1, 1}, {0, 4, 1}, {2, 4, 1}};
    C.n_rows = m;
    C.n_cols = n;
    C.val = {1, 2, 1, 2, 1};
    C.col_ind = {2, 4, 1, 1, 4};
    C.row_ptr = {0, 2, 3, 5};
    X.resize(m, n);
    X << 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1;
  }

  const std::size_t m = 3, n = 5;
  TripletMatrix<double> T;
  CRSMatrix<double> C;
  Eigen::MatrixXd X;
};

TestData data;
constexpr double eps = 1e-10;

TEST_SUITE("TripletToCRS") {
  TEST_CASE("Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> densify" *
            doctest::description("densify() functions")) {
    Eigen::MatrixXd T = densify_TEST(data.T);
    Eigen::MatrixXd C = densify_TEST(data.C);
    REQUIRE(T.rows() == data.T.n_rows);
    REQUIRE(T.cols() == data.T.n_cols);
    REQUIRE(C.rows() == data.C.n_rows);
    REQUIRE(C.cols() == data.C.n_cols);
    CHECK((T - data.X).norm() == doctest::Approx(0.).epsilon(eps));
    CHECK((C - data.X).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("CRSMatrix<SCALAR> tripletToCRS" *
            doctest::description("conversion")) {
    CRSMatrix<double> C = tripletToCRS_TEST(data.T);
    REQUIRE(C.n_rows == data.C.n_rows);
    REQUIRE(C.n_cols == data.C.n_cols);

    CHECK((densify(C) - data.X).norm() == doctest::Approx(0.).epsilon(eps));
  }
}
