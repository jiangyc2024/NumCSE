#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
    n = 8;
    Xstop = Eigen::MatrixXd(2, n);
    Xstop << 1, 3, 3, 1, -1, -3, -3, -1, -3, -1, 1, 3, 3, 1, -1, -3;
    Xpriority = Eigen::MatrixXd(2, n);
    Xpriority << 0, 3, 0, -3, 0, 2.5, 0, -2.5, -3, 0, 3, 0, -2.5, 0, 2.5, 0;

    P1 = Eigen::MatrixXd(2, n);
    P1 << 0.23657, 1.35369, -0.13624, -1.33702, 0.0989619, 0.993235, -0.0735973,
        -1.11657, -2.76114, -2.60103, 2.90403, 2.66831, -2.44302, -2.04656,
        2.31922, 2.20296;
    P2 = Eigen::MatrixXd(2, n);
    P2 << -1.12783, -1.75868, -1.40935, -0.0664574, 1.09654, 1.75873, 1.5195,
        -0.0607661, 1.72169, 0.344036, -0.889686, -1.87847, -1.57535, -0.41511,
        0.834371, 1.88514;

    P3 = Eigen::MatrixXd(2, n);
    P3 << -1.23988, -0.731643, 0.00492048, 1.08039, 1.34128, 0.670982,
        -0.101797, -1.02859, 1.53076, 2.02881, 1.36163, -0.340912, -1.47697,
        -1.99975, -1.47947, 0.374859;
  }

  unsigned int n;
  Eigen::Eigen::MatrixXd Xstop, Xpriority, P1, P2, P3;
};

TestData data;

TEST_SUITE("ShapeIdent") {
  TEST_CASE("Eigen::MatrixXd shape_ident_matrix" *
            doctest::description("shape_ident_matrix()")) {
    Eigen::Eigen::MatrixXd sol, stud;

    sol = shape_ident_matrix(data.Xstop);
    stud = shape_ident_matrix_TEST(data.Xstop);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));
  }
  TEST_CASE("double solve_lsq" * doctest::description("solve_lsq()")) {
    Eigen::MatrixXd A_stud, A_sol;
    double sol, stud;

    sol = solve_lsq(data.Xstop, data.P1, A_sol);
    stud = solve_lsq_TEST(data.Xstop, data.P1, A_stud);

    REQUIRE(A_sol.size() == A_stud.size());
    CHECK((A_sol - A_stud).norm() == doctest::Approx(0.).epsilon(1e-9));
    CHECK(sol - stud == doctest::Approx(0.).epsilon(1e-9));
  }
  TEST_CASE("Shape identify" * doctest::description("identify()")) {
    Eigen::MatrixXd A_stud1, A_sol1, A_stud2, A_sol2, A_stud3, A_sol3;
    Shape s_stud1, s_sol1, s_stud2, s_sol2, s_stud3, s_sol3;

    s_sol1 = identify(data.Xstop, data.Xpriority, data.P1, A_sol1);
    s_stud1 = identify_TEST(data.Xstop, data.Xpriority, data.P1, A_stud1);

    s_sol2 = identify(data.Xstop, data.Xpriority, data.P2, A_sol2);
    s_stud2 = identify_TEST(data.Xstop, data.Xpriority, data.P2, A_stud2);
    s_sol3 = identify(data.Xstop, data.Xpriority, data.P3, A_sol3);
    s_stud3 = identify_TEST(data.Xstop, data.Xpriority, data.P3, A_stud3);

    REQUIRE(A_sol1.size() == A_stud1.size());
    REQUIRE(A_sol2.size() == A_stud2.size());
    REQUIRE(A_sol3.size() == A_stud3.size());
    CHECK((A_sol1 - A_stud1).norm() == doctest::Approx(0.).epsilon(1e-9));
    CHECK((A_sol2 - A_stud2).norm() == doctest::Approx(0.).epsilon(1e-9));
    CHECK((A_sol3 - A_stud3).norm() == doctest::Approx(0.).epsilon(1e-9));
    CHECK(s_sol1 == s_stud1);
    CHECK(s_sol2 == s_stud2);
    CHECK(s_sol3 == s_stud3);
  }
}
