#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

TEST_SUITE("FFTLSQ") {
  TEST_CASE("Eigen::VectorXd eval_p" * doctest::description("helper") *
            doctest::skip()) {}

  TEST_CASE("bool testNormEqMatrix" *
            doctest::description("testNormEqMatrix()")) {
    unsigned int n = 10;
    unsigned int m = 3;
    CHECK(testNormEqMatrix_TEST(n, m));
  }

  TEST_CASE("Eigen::VectorXd find_c" * doctest::description("find_c()")) {
    Eigen::VectorXd d(10);
    d.setRandom();

    auto sol = find_c(d, 3);
    auto stud = find_c_TEST(d, 3);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));

    d << 0.987214, 1.03579, 0.997689, 0.917471, 1.00474, 0.92209, 1.03517,
        1.08863, 0.904992, 0.956089;

    sol = find_c(d, 3);
    stud = find_c_TEST(d, 3);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  // clang-format off
  TEST_CASE("void fitEllipse" *
      doctest::description("fit polynomial coefficients and plot ellipses")) {
    // clang-format on
    MESSAGE("This function wasn't tested. Run the program and check the plot.");
  }
}
