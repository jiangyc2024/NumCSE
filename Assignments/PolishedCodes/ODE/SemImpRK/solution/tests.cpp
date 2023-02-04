#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

typedef Eigen::Matrix<double, 1, 1> Vec;

struct TestData {
  TestData() {
    T = 1.1;
    n = 30;
    Y0.resize(1, 1);
    Y0 << 1;
  }

  double T;
  unsigned int n;
  Eigen::MatrixXd Y0;
};

TestData data;

Vec f(Vec x) { return x; }

Vec df(Vec x) {
  Vec out;
  out << 1;
  return out;
}

TEST_SUITE("SemImpRK") {
  TEST_CASE("std::vector<Eigen::VectorXd> solveRosenbrock" *
            doctest::description("Check solution at final time")) {
    auto sol = solveRosenbrock(f, df, data.Y0, data.n, data.T);
    auto stud = solveRosenbrock_TEST(f, df, data.Y0, data.n, data.T);

    const bool samesize = sol.back().rows() == stud.back().rows() &&
                          sol.back().cols() == stud.back().cols();
    REQUIRE(samesize);
    CHECK((sol.back() - stud.back()).norm() ==
          doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("double cvgRosenbrock" * doctest::description("cvgRosenbrock()")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
