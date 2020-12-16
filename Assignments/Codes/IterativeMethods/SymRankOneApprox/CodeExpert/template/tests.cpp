#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

using namespace Eigen;

struct TestData {
  TestData() {
    eps = 1e-4;
    M << 1 + eps, 2, 1, 2, 1 + eps, 1, 1, 1, 1 + eps;

    v << 2, 0, 1;
    b << 1, 2, 3, 2, 3, 4, 3, 4, 5;

    MM << 1, 0, 1, 1, 0, 0, 2, 1, 0;
    rtol = 1e-6;
    atol = 1e-8;

    z_bestApprox << 1.21315, 1.21315, 0.888086;
    w_kronProd << 10, 16, 12;
    z_approx << 1.16657, 0.441537, 0.859084;
  }

  // input data for approximation
  // symRankOneBestApproxSym
  MatrixXd M;
  double eps;

  // computeKronProdVecMult
  VectorXd v;
  VectorXd b;

  // symmRankOneApprox
  MatrixXd MM;
  double rtol;
  double atol;

  // exact solution to test against
  VectorXd z_bestApprox;
  VectorXd w_kronProd;
  VectorXd z_approx;
};

TestData data;

TEST_SUITE("SymRankOneApprox") {
  TEST_CASE("VectorXd symRankOneBestApproxSym" * doctest::description("approximation of Forbenius norm for a symmetric matrixW")) {
    VectorXd sol = data.z_bestApprox;
    VectorXd stud = symRankOneBestApproxSym_TEST(data.M);

    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-4));
  }

  TEST_CASE("VectorXd computeKronProdVecMult" * doctest::description("evaluation of the Kronecker product")) {
    VectorXd sol = data.w_kronProd;
    VectorXd stud = computeKronProdVecMult_TEST(data.v, data.b);

    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-4));
  }

  TEST_CASE("VectorXd symmRankOneApprox" * doctest::description("approximation of Forbenius norm using the Gauss-Newton iteration")) {
    VectorXd sol = data.z_approx;
    VectorXd stud = symmRankOneApprox_TEST(data.MM, data.rtol, data.atol);

    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-4));
  }
}
