#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "compactstorageqr.hpp"
#include "solution.hpp"
#include <utility>
#include <cmath>

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
  }
};

TestData data;

TEST_SUITE("CompactStorageQR") {

  TEST_CASE("matmult()" * doctest::description("Testing output, not efficiency")) {
    // Comparing result of matmult with solution
    // Test case 1
    unsigned N = 7;
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(N,N);
    M /= 10;
    CompactStorageQR obj(M.data(),N);
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(N,N);
    Eigen::MatrixXd res = obj.matmult(X);

    solution::CompactStorageQR objsol(M.data(),N);
    Eigen::MatrixXd ressol = objsol.matmult(X);

    CHECK( (res-ressol).norm() == doctest::Approx(0).epsilon(1e-10) );

    // Test case 2
    unsigned n = 9;
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(n,n);
    m /= 10;
    CompactStorageQR obj1(m.data(),n);
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(n,n);
    Eigen::MatrixXd res1 = obj1.matmult(x);

    solution::CompactStorageQR objsol1(m.data(),n);
    Eigen::MatrixXd ressol1 = objsol1.matmult(x);

    CHECK( (res1-ressol1).norm() == doctest::Approx(0).epsilon(1e-10) );
  }

  TEST_CASE("solve()" * doctest::description("Testing output, not efficiency")) {
    // Comparing result of solve with solution
    // Test case 1
    unsigned N = 7;
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(N,N);
    M /= 10;
    CompactStorageQR obj(M.data(),N);
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(N,N);
    Eigen::MatrixXd res = obj.solve(X);

    solution::CompactStorageQR objsol(M.data(),N);
    Eigen::MatrixXd ressol = objsol.solve(X);

    CHECK( (res-ressol).norm() == doctest::Approx(0).epsilon(1e-10) );

    // Test case 2
    unsigned n = 9;
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(n,n);
    m /= 10;
    CompactStorageQR obj1(m.data(),n);
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(n,n);
    Eigen::MatrixXd res1 = obj1.solve(x);

    solution::CompactStorageQR objsol1(m.data(),n);
    Eigen::MatrixXd ressol1 = objsol1.solve(x);

    CHECK( (res1-ressol1).norm() == doctest::Approx(0).epsilon(1e-10) );
  }

  TEST_CASE("det()" * doctest::description("Testing output")) {
    unsigned N = 7;
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(N,N);
    M /= 10;
    CompactStorageQR obj(M.data(),N);
    double res = obj.det();

    solution::CompactStorageQR objsol(M.data(),N);
    double ressol = objsol.det();
    CHECK( std::abs(res-ressol)/std::abs(ressol) == doctest::Approx(0).epsilon(1e-10) );

    unsigned n = 9;
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(n,n);
    m /= 10;
    CompactStorageQR obj1(m.data(),n);
    double res1 = obj1.det();

    solution::CompactStorageQR objsol1(m.data(),n);
    double ressol1 = objsol1.det();

    CHECK( std::abs(res1-ressol1)/std::abs(ressol1) == doctest::Approx(0).epsilon(1e-10) );
  }
}
