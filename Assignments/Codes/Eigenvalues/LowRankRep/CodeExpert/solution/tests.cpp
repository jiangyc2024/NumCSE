#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
    m = 3;
    n = 2;
    k = 2;
    X.resize(m, n);
    X << 5, 0, 2, 1, 7, 4;
    A.resize(m, k);
    B.resize(n, k);
    A << 2, 1, 2, 3, 6, 1;
    B << 4, 4, 5, 0;
    Ax.resize(m, k);
    Ay.resize(m, k);
    Bx.resize(n, k);
    By.resize(n, k);
    Ax << 1, 0, 9, 2, 6, 3;
    Ay << 8, -2, 3, 4, 5, 8;
    Bx << 2, 1, 2, 3;
    By << 4, 4, 5, 0;
  }

  std::size_t m, n, k;
  Eigen::MatrixXd X, A, B, Ax, Ay, Bx, By;
};

TestData data;

TEST_SUITE("LowRankRep") {
  TEST_CASE("std::pair<MatrixXd,MatrixXd> factorize_X_AB" *
            doctest::description("factor")) {
    auto sol = factorize_X_AB(data.X, data.k);
    auto stud = factorize_X_AB_TEST(data.X, data.k);

    REQUIRE(sol.first.rows() == stud.first.rows());
    REQUIRE(sol.second.rows() == stud.second.rows());
    REQUIRE(sol.first.cols() == stud.first.cols());
    REQUIRE(sol.second.cols() == stud.second.cols());
    CHECK((sol.first - stud.first).norm() == doctest::Approx(0.).epsilon(1e-6));
    CHECK((sol.second - stud.second).norm() ==
          doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("std::tuple<MatrixXd, MatrixXd, MatrixXd> svd_AB" *
            doctest::description("svd U, S, V")) {
    auto sol = svd_AB(data.A, data.B);
    auto stud = svd_AB_TEST(data.A, data.B);

    auto sol0 = std::get<0>(sol);
    auto stud0 = std::get<0>(stud);
    auto sol1 = std::get<1>(sol);
    auto stud1 = std::get<1>(stud);
    auto sol2 = std::get<2>(sol);
    auto stud2 = std::get<2>(stud) REQUIRE(sol0.rows() == stud0.rows());
    REQUIRE(sol0.cols() == stud0.cols());
    REQUIRE(sol1.rows() == stud1.rows());
    REQUIRE(sol1.cols() == stud1.cols());
    REQUIRE(sol2.rows() == stud2.rows());
    REQUIRE(sol2.cols() == stud2.cols());
    CHECK((sol0 - stud0).norm() == doctest::Approx(0.).epsilon(1e-6));
    CHECK((sol1 - stud1).norm() == doctest::Approx(0.).epsilon(1e-6));
    CHECK((sol2 - stud2).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("std::pair<MatrixXd,MatrixXd> rank_k_approx" *
            doctest::description("approx Z")) {
    auto sol = rank_k_approx(data.Ax, data.Ay, data.Bx, data.By);
    auto stud = rank_k_approx_TEST(data.Ax, data.Ay, data.Bx, data.By);

    REQUIRE(sol.first.rows() == stud.first.rows());
    REQUIRE(sol.second.rows() == stud.second.rows());
    REQUIRE(sol.first.cols() == stud.first.cols());
    REQUIRE(sol.second.cols() == stud.second.cols());
    CHECK((sol.first - stud.first).norm() == doctest::Approx(0.).epsilon(1e-6));
    CHECK((sol.second - stud.second).norm() ==
          doctest::Approx(0.).epsilon(1e-6));
  }
}
