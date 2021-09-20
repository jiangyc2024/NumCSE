#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() { A = Eigen::Matrix3d::Identity(); }

  Eigen::Matrix3d A;
};

TestData data;
constexpr double eps = 1e-10;

TEST_SUITE("MatrixReduce") {
  TEST_CASE("double average" * doctest::description("average")) {
    const double sol = average(data.A);
    const double stud = average_TEST(data.A);

    CHECK(sol == doctest::Approx(stud).epsilon(eps));
  }

  TEST_CASE("double percent_zero" * doctest::description("percent_zero")) {
    const double sol = percent_zero(data.A);
    const double stud = percent_zero_TEST(data.A);

    CHECK(sol == doctest::Approx(stud).epsilon(eps));
  }

  TEST_CASE("bool has_zero_column" * doctest::description("has_zero_column")) {
    Eigen::MatrixXd B = Eigen::MatrixXd::Identity(4, 5);

    CHECK(has_zero_column(B) == has_zero_column_TEST(B));
    CHECK(has_zero_column(B.transpose()) ==
          has_zero_column_TEST(B.transpose()));
  }

  TEST_CASE("Eigen::MatrixXd columns_sum_to_zero" *
            doctest::description("columns_sum_to_zero")) {
    Eigen::Matrix3d C =
        Eigen::Matrix3d::Random() + Eigen::Matrix3d::Constant(1);

    auto sol = columns_sum_to_zero(C);
    auto stud = columns_sum_to_zero_TEST(C);
    REQUIRE(sol.rows() == stud.rows());
    REQUIRE(sol.cols() == stud.cols());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(eps));
  }
}
