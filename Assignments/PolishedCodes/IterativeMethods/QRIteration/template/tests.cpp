#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

// includes for test data
#include "qriteration.hpp"
#include <Eigen/Dense>
#include <cmath>

TEST_SUITE("QR-Iteration") {

  TEST_CASE("qrStep()" *
            doctest::description(
                "Checking the returned vectors for a random symmetric "
                "tridiagonal matrix (Not testing efficiency)")) {
    unsigned n = 5;
    Eigen::VectorXd diag = Eigen::VectorXd::Random(n);
    Eigen::VectorXd sub_diag = Eigen::VectorXd::Random(n - 1);

    Eigen::MatrixXd TriDiagMat = Eigen::MatrixXd::Constant(n, n, 0);
    TriDiagMat.diagonal() = diag;
    TriDiagMat.diagonal(1) = sub_diag;
    TriDiagMat.diagonal(-1) = sub_diag;

    Eigen::HouseholderQR<Eigen::MatrixXd> qr(TriDiagMat);
    Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd T_p = Q.transpose() * TriDiagMat * Q;

    Eigen::VectorXd diag_p = T_p.diagonal();
    Eigen::VectorXd sub_diag_p = T_p.diagonal(1);

    Eigen::VectorXd v1, v2;
    std::tie(v1, v2) = qrStep(diag, sub_diag);

    CHECK((v1.cwiseAbs() - diag_p.cwiseAbs()).norm() ==
          doctest::Approx(0.).epsilon(1e-10));
    CHECK((v2.cwiseAbs() - sub_diag_p.cwiseAbs()).norm() ==
          doctest::Approx(0.).epsilon(1e-10));
  }
}