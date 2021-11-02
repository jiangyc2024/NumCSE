#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
	TestData() {
		A.resize(5, 5);
		A << 10, 2, 3, 4, 5, 6, 20, 8, 9, 1, 1, 2, 30, 4, 5, 6, 7, 8, 20, 0, 1, 2, 3, 4, 10;
	}
	
	Eigen::MatrixXd A;
};

TestData data;
constexpr double eps = 1e-12;

TEST_SUITE("Lyapunov") {
	TEST_CASE("Eigen::SparseMatrix<double> buildC" * doctest::description("System matrix")) {
		auto C_sol = buildC(data.A);
		auto C_stud = buildC_TEST(data.A);
		
		REQUIRE(C_sol.rows() == C_stud.rows());
		REQUIRE(C_sol.cols() == C_stud.cols());
		CHECK((C_sol - C_stud).norm() == doctest::Approx(0.).epsilon(eps));
	}
	
	TEST_CASE("void solveLyapunov" * doctest::description("Solution")) {
		Eigen::MatrixXd X_sol, X_stud;
		
		solveLyapunov(data.A, X_sol);
		solveLyapunov_TEST(data.A, X_stud);
		
		REQUIRE(X_sol.rows() == X_stud.rows());
		REQUIRE(X_sol.cols() == X_stud.cols());
		CHECK((X_sol - X_stud).norm() == doctest::Approx(0.).epsilon(eps));
	}
	
	TEST_CASE("bool testLyapunov" * doctest::description("test")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
	
}

