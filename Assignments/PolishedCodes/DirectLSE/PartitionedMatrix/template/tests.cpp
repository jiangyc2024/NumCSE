#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
    const unsigned int n = 25;
    std::srand(n);
    R = Eigen::MatrixXd::Random(n, n).triangularView<Eigen::Upper>();
    v = Eigen::VectorXd::Random(n);
    u = Eigen::VectorXd::Random(n);
    bb = Eigen::VectorXd::Random(n + 1);
  }

  Eigen::MatrixXd R;
  Eigen::VectorXd v, u, bb;
};

TestData data;
constexpr double eps = 1e-12;

TEST_SUITE("PartitionedMatrix") {
  TEST_CASE("void solvelse" * doctest::description("solvelse()")) {
    Eigen::VectorXd xo_sol, xo_stud;

    solvelse(data.R, data.v, data.u, data.bb, xo_sol);
    solvelse_TEST(data.R, data.v, data.u, data.bb, xo_stud);

    REQUIRE(xo_sol.size() == xo_stud.size());
    CHECK((xo_sol - xo_stud).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("bool testSolveLSE" * doctest::description("testSolveLSE()")) {
    Eigen::VectorXd xe_sol, xe_stud;

    bool check_sol = testSolveLSE(data.R, data.v, data.u, data.bb, xe_sol);
    bool check_stud =
        testSolveLSE_TEST(data.R, data.v, data.u, data.bb, xe_stud);

    // check if LU-Decomposition was done correctly
    REQUIRE(xe_sol.size() == xe_stud.size());
    CHECK((xe_sol - xe_stud).norm() == doctest::Approx(0.).epsilon(eps));
    CHECK(check_stud == check_sol);
  }
}
