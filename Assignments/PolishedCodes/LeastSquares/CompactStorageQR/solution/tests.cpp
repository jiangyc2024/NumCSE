#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <utility>

#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

TEST_SUITE("CompactStorageQR") {
  TEST_CASE("Eigen::MatrixXd matmult [OUT OF CLASS]" *
            doctest::description("Testing output, not efficiency")) {
    // Comparing result of matmult with solution
    // Test case 1
    constexpr unsigned int N = 7;
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(N, N);
    M /= 10;
    CompactStorageQR_TEST obj(M.data(), N);
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(N, N);
    Eigen::MatrixXd res = obj.matmult(X);

    CompactStorageQR objsol(M.data(), N);
    Eigen::MatrixXd ressol = objsol.matmult(X);

    REQUIRE(res.rows() == ressol.rows());
    REQUIRE(res.cols() == ressol.cols());
    CHECK((res - ressol).norm() == doctest::Approx(0).epsilon(1e-10));

    // Test case 2
    constexpr unsigned int n = 9;
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(n, n);
    m /= 10;
    CompactStorageQR_TEST obj1(m.data(), n);
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd res1 = obj1.matmult(x);

    CompactStorageQR objsol1(m.data(), n);
    Eigen::MatrixXd ressol1 = objsol1.matmult(x);

    REQUIRE(res1.rows() == ressol1.rows());
    REQUIRE(res1.cols() == ressol1.cols());
    CHECK((res1 - ressol1).norm() == doctest::Approx(0).epsilon(1e-10));
  }

  TEST_CASE("Eigen::MatrixXd solve [OUT OF CLASS]" *
            doctest::description("Testing output, not efficiency")) {
    // Comparing result of solve with solution
    // Test case 1
    constexpr unsigned int N = 7;
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(N, N);
    M /= 10;
    CompactStorageQR_TEST obj(M.data(), N);
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(N, N);
    Eigen::MatrixXd res = obj.solve(X);

    CompactStorageQR objsol(M.data(), N);
    Eigen::MatrixXd ressol = objsol.solve(X);

    REQUIRE(res.rows() == ressol.rows());
    REQUIRE(res.cols() == ressol.cols());
    CHECK((res - ressol).norm() == doctest::Approx(0).epsilon(1e-10));

    // Test case 2
    constexpr unsigned int n = 9;
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(n, n);
    m /= 10;
    CompactStorageQR_TEST obj1(m.data(), n);
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd res1 = obj1.solve(x);

    CompactStorageQR objsol1(m.data(), n);
    Eigen::MatrixXd ressol1 = objsol1.solve(x);

    REQUIRE(res1.rows() == ressol1.rows());
    REQUIRE(res1.cols() == ressol1.cols());
    CHECK((res1 - ressol1).norm() == doctest::Approx(0).epsilon(1e-10));
  }

  TEST_CASE("double det [OUT OF CLASS]" *
            doctest::description("Testing output")) {
    constexpr unsigned int N = 7;
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(N, N);
    M /= 10;
    CompactStorageQR_TEST obj(M.data(), N);
    const double res = obj.det();

    CompactStorageQR objsol(M.data(), N);
    const double ressol = objsol.det();
    CHECK(std::abs(res - ressol) / std::abs(ressol) ==
          doctest::Approx(0).epsilon(1e-10));

    constexpr unsigned int n = 9;
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(n, n);
    m /= 10;
    CompactStorageQR_TEST obj1(m.data(), n);
    const double res1 = obj1.det();

    CompactStorageQR objsol1(m.data(), n);
    const double ressol1 = objsol1.det();

    CHECK(std::abs(res1 - ressol1) / std::abs(ressol1) ==
          doctest::Approx(0).epsilon(1e-10));
  }
}
