#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Core>
#include <Eigen/SparseCore>

// We use this as type of each index (will be int or similar)
using index_t = Eigen::Index;
// Use this to contain the "size" of our matrix
using shape_t = std::array<index_t, 2>;

struct TestData {
  TestData() {
    n = 16;
    m = 16;
    S << 0, 1, 0, 1, -4, 1, 0, 1, 0;
    f = define_f(n, m);
  }

  index_t n, m;
  Eigen::Matrix3d S;
  std::function<double(index_t, index_t)> f;
};

TestData data;

TEST_SUITE("GridFunction") {
  TEST_CASE("std::function<double(index_t, index_t)> define_f" *
            doctest::description("skipped") * doctest::skip()) {}

  TEST_CASE("void eval" * doctest::description("test matrix")) {
    Eigen::MatrixXd X_sol(data.n, data.m);
    Eigen::MatrixXd X_stud(data.n, data.m);
    eval(X_sol, data.f);
    eval_TEST(X_stud, data.f);

    REQUIRE(X_sol.rows() == X_stud.rows());
    REQUIRE(X_sol.cols() == X_stud.cols());
    CHECK((X_sol - X_stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("inline index_t to_vector_index" * doctest::description("skipped") *
            doctest::skip()) {}

  TEST_CASE("Eigen::SparseMatrix<double> build_matrix" *
            doctest::description("build_matrix")) {
    Eigen::SparseMatrix<double> A_small_sol =
        build_matrix(data.S, shape_t{3, 3});
    Eigen::SparseMatrix<double> A_small_stud =
        build_matrix_TEST(data.S, shape_t{3, 3});

    REQUIRE(A_small_sol.rows() == A_small_stud.rows());
    REQUIRE(A_small_sol.cols() == A_small_stud.cols());
    CHECK((A_small_stud - A_small_sol).norm() ==
          doctest::Approx(0.).epsilon(1e-8));

    Eigen::SparseMatrix<double> A_sol =
        build_matrix(data.S, shape_t{data.n, data.m});
    Eigen::SparseMatrix<double> A_stud =
        build_matrix_TEST(data.S, shape_t{data.n, data.m});

    REQUIRE(A_sol.rows() == A_stud.rows());
    REQUIRE(A_sol.cols() == A_stud.cols());
    CHECK((A_stud - A_sol).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("void mult" * doctest::description("mult")) {
    Eigen::SparseMatrix<double> A =
        build_matrix(data.S, shape_t{data.n, data.m});
    Eigen::MatrixXd X(data.n, data.m);
    eval(X, data.f);
    Eigen::MatrixXd Y_sol, Y_stud;

    mult(A, X, Y_sol);
    mult_TEST(A, X, Y_stud);

    REQUIRE(Y_sol.rows() == Y_stud.rows());
    REQUIRE(Y_sol.cols() == Y_stud.cols());
    CHECK((Y_sol - Y_stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("void solve" * doctest::description("solve")) {
    Eigen::SparseMatrix<double> A =
        build_matrix(data.S, shape_t{data.n, data.m});
    Eigen::MatrixXd X(data.n, data.m);
    eval(X, data.f);
    Eigen::MatrixXd Y;
    mult(A, X, Y);
    Eigen::MatrixXd X_sol, X_stud;

    solve(A, Y, X_sol);
    solve_TEST(A, Y, X_stud);

    REQUIRE(X_sol.rows() == X_stud.rows());
    REQUIRE(X_sol.cols() == X_stud.cols());
    CHECK((X_sol - X_stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }
}
