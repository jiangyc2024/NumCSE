#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
    X.resize(3, 2);
    X << 5, 0, 2, 1, 7, 4;
    A.resize(m, k);
    B.resize(n, k);
    A << 2, 1, 2, 3, 6, 1;
    B << 4, 4, 5, 0, 7, 4, 3, 2;
    Ax.resize(6, k);
    Ay.resize(6, k);
    Bx.resize(5, k);
    By.resize(5, k);
    Ax << 0, 9, 2, 6, 3, 5, -2, 3, 4, 8, 9, 0;
    Ay << 8, -2, 3, 4, 5, 8, 7, 5, 6, 3, 2, 1;
    Bx << 2, 1, 2, -3, 6, 7, 8, 4, 5, 1;
    By << 4, 4, -5, 0, 3, 2, 5, 9, 0, 5;
  }

  const std::size_t k = 2;
  const std::size_t m = 3;
  const std::size_t n = 4;
  Eigen::MatrixXd X, A, B, Ax, Ay, Bx, By;
};

TestData data;
constexpr double eps = 1e-6;

TEST_SUITE("LowRankRep") {
  TEST_CASE("std::pair<Eigen::MatrixXd, Eigen::MatrixXd> factorize_X_AB" *
            doctest::description("factor")) {
    auto stud = factorize_X_AB_TEST(data.X, data.k);

    REQUIRE(stud.first.rows() == 3);
    REQUIRE(stud.first.cols() == data.k);
    REQUIRE(stud.second.rows() == 2);
    REQUIRE(stud.second.cols() == data.k);
    Eigen::MatrixXd X = stud.first * stud.second.transpose();
    CHECK((X - data.X).norm() == doctest::Approx(0.).epsilon(eps));
  }

  // clang-format off
  TEST_CASE("std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> svd_AB" *
      doctest::description("svd U, S, V")) {
    // clang-format on
    auto sol = svd_AB(data.A, data.B);
    auto stud = svd_AB_TEST(data.A, data.B);

    auto sol0 = std::get<0>(sol);
    auto stud0 = std::get<0>(stud);
    auto sol1 = std::get<1>(sol);
    auto stud1 = std::get<1>(stud);
    auto sol2 = std::get<2>(sol);
    auto stud2 = std::get<2>(stud);
    REQUIRE(sol0.rows() == stud0.rows());
    REQUIRE(sol0.cols() == stud0.cols());
    REQUIRE(sol1.rows() == stud1.rows());
    REQUIRE(sol1.cols() == stud1.cols());
    REQUIRE(sol2.rows() == stud2.rows());
    REQUIRE(sol2.cols() == stud2.cols());
    CHECK((sol0 - stud0).norm() == doctest::Approx(0.).epsilon(eps));
    CHECK((sol1 - stud1).norm() == doctest::Approx(0.).epsilon(eps));
    CHECK((sol2 - stud2).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("std::pair<Eigen::MatrixXd, Eigen::MatrixXd> rank_k_approx" *
            doctest::description("approx Z")) {
    auto sol = rank_k_approx(data.Ax, data.Ay, data.Bx, data.By);
    auto stud = rank_k_approx_TEST(data.Ax, data.Ay, data.Bx, data.By);

    REQUIRE(sol.first.rows() == stud.first.rows());
    REQUIRE(sol.second.rows() == stud.second.rows());
    REQUIRE(sol.first.cols() == stud.first.cols());
    REQUIRE(sol.second.cols() == stud.second.cols());
    CHECK((sol.first - stud.first).norm() == doctest::Approx(0.).epsilon(eps));
    CHECK((sol.second - stud.second).norm() ==
          doctest::Approx(0.).epsilon(eps));
  }
}
