#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    T = 1;
    n = 3;
    Y0.resize(n, n);
    Y0 << 1, 1, 0, 0, 3, 2, 1, 5, 2;
  }

  double T;
  unsigned int n;
  Eigen::MatrixXd Y0;
};

TestData data;

TEST_SUITE("NLMatODE") {
  TEST_CASE("Eigen::MatrixXd matode" *
            doctest::description("Check matode matrix at T")) {
    const Eigen::MatrixXd sol = matode(data.Y0, data.T);
    const Eigen::MatrixXd stud = matode_TEST(data.Y0, data.T);

    const bool samesize =
        sol.rows() == stud.rows() && sol.cols() == stud.cols();
    REQUIRE(samesize);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("bool checkinvariant" *
            doctest::description(
                "Test whether invariant for Y0 was preserved or not")) {
    const bool sol = checkinvariant(data.Y0, data.T);
    const bool stud = checkinvariant_TEST(data.Y0, data.T);

    CHECK(sol == stud);
  }

  TEST_CASE("double cvgDiscreteGradientMethod" *
            doctest::description("Test fitted convergance rate")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
