#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    A = Eigen::MatrixXd::Random(10, 10);
    B = Eigen::MatrixXd::Random(10, 10);
  }
  Eigen::MatrixXd A;
  Eigen::MatrixXd B;
};

TestData data;

TEST_SUITE("Gram-Schmidt") {
  TEST_CASE("Eigen::MatrixXd gram_schmidt" *
            doctest::description("gram_schmidt(A)")) {
    Eigen::MatrixXd C_sol, D_sol;
    Eigen::MatrixXd C_stud, D_stud;
    C_sol = gram_schmidt(data.A);
    C_stud = gram_schmidt_TEST(data.A);
    D_sol = gram_schmidt(data.B);
    D_stud = gram_schmidt_TEST(data.B);

    REQUIRE(D_sol.rows() == D_stud.rows());
    REQUIRE(D_sol.cols() == D_stud.cols());
    REQUIRE(C_sol.rows() == C_stud.rows());
    REQUIRE(C_sol.cols() == C_stud.cols());

    CHECK((C_sol - C_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
    CHECK((D_sol - D_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("bool testGramSchmidt" * doctest::description("Test by student")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
