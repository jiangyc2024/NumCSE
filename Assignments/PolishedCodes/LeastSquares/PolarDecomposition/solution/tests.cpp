#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>
#include <cmath>

#include "polardecomposition.hpp"

TEST_SUITE("Polar Decomposition") {
  TEST_CASE("void initialize" *
            doctest::description("Checking the initialized factors Q and M")) {
    constexpr unsigned int m = 5;
    constexpr unsigned int n = 4;
    // Random m X n matrix
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(m, n);
    PolarDecomposition polardecomp(X);
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd X_recreated = Id;
    polardecomp.applyM(X_recreated);
    // Checking if the computed M factor is symmetric
    CHECK((X_recreated - X_recreated.transpose()).norm() ==
          doctest::Approx(0.).epsilon(1e-10));
    polardecomp.applyQ(X_recreated);
    // Checking if X == QM
    CHECK((X - X_recreated).norm() == doctest::Approx(0.).epsilon(1e-10));
    Eigen::MatrixXd Q = Id;
    polardecomp.applyQ(Q);
    // Checking if Q has orthogonal columns
    CHECK((Q.transpose() * Q - Id).norm() ==
          doctest::Approx(0.).epsilon(1e-10));
  }

  TEST_CASE("PolarDecomposition" *
            doctest::description(
                "Checking the factors Q and M, not the efficiency")) {
    constexpr unsigned int m = 11;
    constexpr unsigned int n = 7;
    constexpr unsigned int k = 3;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(m, k);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(n, k);
    Eigen::MatrixXd X = A * B.transpose();
    PolarDecomposition polardecomp(A, B);
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd X_recreated = Id;
    polardecomp.applyM(X_recreated);
    // Checking symmetry of the factor M
    CHECK((X_recreated - X_recreated.transpose()).norm() ==
          doctest::Approx(0.).epsilon(1e-10));
    polardecomp.applyQ(X_recreated);
    // Checking if X = QM
    CHECK((X - X_recreated).norm() == doctest::Approx(0.).epsilon(1e-10));
    Eigen::MatrixXd Q = Id;
    polardecomp.applyQ(Q);
    // Checking Q'Q = Id
    CHECK((Q.transpose() * Q - Id).norm() ==
          doctest::Approx(0.).epsilon(1e-10));
  }
}
