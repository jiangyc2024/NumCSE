#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "doctest.h"
#include "matPow.hpp"

TEST_SUITE("MatPow") {
  TEST_CASE("Eigen::MatrixXcd matPow" * doctest::description("matrix power")) {
    Eigen::MatrixXcd A = construct_matrix(20);
    Eigen::MatrixXcd sol = A.pow(5);
    Eigen::MatrixXcd stud = matPow(A, 5);
    const bool checksize =
        sol.rows() == stud.rows() && sol.cols() == stud.cols();
    REQUIRE(checksize);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::MatrixXcd construct_matrix" *
            doctest::description("matrix construction") * doctest::skip()) {}

  TEST_CASE("void tabulateRuntime" *
            doctest::description("Runtime measurements")) {
    MESSAGE("This function isn't tested. Run the program to see its output.");
  }
}