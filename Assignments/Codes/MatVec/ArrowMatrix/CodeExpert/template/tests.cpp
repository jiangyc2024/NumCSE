#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"
#include <Eigen/Dense>

TEST_SUITE("Arrow matrix x vector multiplication") {
	
	TEST_CASE("void arrow_matrix_2_times_x" * doctest::description("already implemented") * doctest::skip()) {}
	
	TEST_CASE("void efficient_arrow_matrix_2_times_x" * doctest::description("efficient arrow matrix x vector multiplication")) {
		
		Eigen::VectorXd a = Eigen::VectorXd::LinSpaced(100, 0, 100);
		Eigen::VectorXd d = Eigen::VectorXd::LinSpaced(100, 100, 0);
		Eigen::VectorXd y_sol, y_stud;
		
		efficient_arrow_matrix_2_times_x(d, a, a, y_sol);
		efficient_arrow_matrix_2_times_x_TEST(d, a, a, y_stud);
		
		REQUIRE(y_sol.size() == y_stud.size());
		CHECK((y_sol - y_stud).norm() == doctest::Approx(0.).epsilon(1e-8));
		
		// test identity matrix
		Eigen::VectorXd ones_d = Eigen::VectorXd::Ones(100);
		Eigen::VectorXd zero_a = Eigen::VectorXd:Zero(100);
		Eigen::VectorXd y_stud_ident;
		
		efficient_arrow_matrix_2_times_x_TEST(ones_d, zero_a, a, y_stud_ident);
		
		REQUIRE(y_stud_ident.size() == a.size());
		CHECK((y_stud_ident - a).norm() == doctest::Approx(0.).epsilon(1e-8));
		
		Eigen::VectorXd a2 = Eigen::VectorXd::Random(200);
		Eigen::VectorXd d2 = Eigen::VectorXd::Random(200);
		Eigen::VectorXd x2 = Eigen::VectorXd::Random(200);
		Eigen::VectorXd y_sol2, y_stud2;
		
		efficient_arrow_matrix_2_times_x(d2, a2, x2, y_sol2);
		efficient_arrow_matrix_2_times_x_TEST(d2, a2, x2, y_stud2);
		
		REQUIRE(y_sol2.size() == y_stud2.size());
		CHECK((y_sol2 - y_stud2).norm() == doctest::Approx(0.).epsilon(1e-8));
		
	}
	
	TEST_CASE("void runtime_arrow_matrix" * doctest::description("Test and compare runtime")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
	
	TEST_CASE("duration_t timing" * doctest::description("implemented only in solution") * doctest::skip()) {}
	
	TEST_CASE("void runtime_arrow_matrix_with_chrono" * doctest::description("implemented only in solution") * doctest::skip()) {}
	
}
