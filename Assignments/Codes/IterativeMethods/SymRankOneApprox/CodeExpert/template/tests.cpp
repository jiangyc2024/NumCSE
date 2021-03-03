#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

using namespace Eigen;

TEST_SUITE("SymRankOneApprox") {
  TEST_CASE("VectorXd symRankOneBestApproxSym" *
            doctest::description(
                "approximation of Frobenius norm for a symmetric matrix")) {
    const double eps = 1e-4;
    MatrixXd M(3, 3);
    M << 1 + eps, 2, 1, 2, 1 + eps, 1, 1, 1, 1 + eps;
    VectorXd z_bestApprox(3);
    z_bestApprox << 1.21315, 1.21315, 0.888086;

    VectorXd stud = symRankOneBestApproxSym_TEST(M);

    REQUIRE(z_bestApprox.size() == stud.size());
    CHECK((z_bestApprox - stud).norm() == doctest::Approx(0.).epsilon(1e-4));
  }

  TEST_CASE("VectorXd computeKronProdVecMult" *
            doctest::description("evaluation of the Kronecker product")) {
    VectorXd v(3);
    v << 2, 0, 1;
    VectorXd b(9);
    b << 1, 2, 3, 2, 3, 4, 3, 4, 5;
    VectorXd w_kronProd(3);
    w_kronProd << 10, 16, 22;

    VectorXd stud = computeKronProdVecMult_TEST(v, b);

    REQUIRE(w_kronProd.size() == stud.size());
    CHECK((w_kronProd - stud).norm() == doctest::Approx(0.).epsilon(1e-4));
  }

  TEST_CASE("VectorXd symmRankOneApprox" *
            doctest::description(
            "approximation of Frobenius norm using the Gauss-Newton iteration")) {
    MatrixXd MM(3, 3);
    MM << 1, 0, 1, 1, 0, 0, 2, 1, 0;
    constexpr double rtol = 1e-6;
    constexpr double atol = 1e-8;
    VectorXd z_approx(3);
    z_approx << 1.16657, 0.441537, 0.859084;

    VectorXd stud = symmRankOneApprox_TEST(MM, rtol, atol);

    REQUIRE(z_approx.size() == stud.size());
    CHECK((z_approx - stud).norm() == doctest::Approx(0.).epsilon(1e-4));
  }
}
