#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

TEST_SUITE("MatrixBlocks") {
  TEST_CASE("MatrixXd zero_row_col" * doctest::description("zero_row_col")) {
    const Eigen::Matrix3d A = Eigen::Matrix3d::Constant(-1);
    const auto sol = zero_row_col(A, 0, 1);
    const auto stud = zero_row_col_TEST(A, 0, 1);

    REQUIRE(sol.rows() == stud.rows());
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("MatrixXd swap_left_right_blocks" *
            doctest::description("swap_left_right_blocks")) {
    const Eigen::MatrixXd B = Eigen::MatrixXd::Identity(4, 3);
    const auto sol = swap_left_right_blocks(B, 2);
    const auto stud = swap_left_right_blocks_TEST(B, 2);

    REQUIRE(sol.rows() == stud.rows());
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("MatrixXd tridiagonal" * doctest::description("tridiagonal")) {
    const auto sol = tridiagonal(4, -1, 2, -1);
    const auto stud = tridiagonal_TEST(4, -1, 2, -1);

    REQUIRE(sol.rows() == stud.rows());
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }
}
