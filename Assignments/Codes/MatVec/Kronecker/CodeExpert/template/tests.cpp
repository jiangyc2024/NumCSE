#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"
#include <Eigen/Dense>

TEST_CASE("Kronecker test") {
	Eigen::MatrixXd A(2, 2);
	Eigen::MatrixXd B(2, 2);
	A << 1, 2, 3, 4;
	B << 5, 6, 7, 8;
	srand(5);
	Eigen::VectorXd x = Eigen::VectorXd::Random(4);
	
	SUBCASE("kron(A,B)*x") {
		Eigen::MatrixXd C_sol;
		Eigen::MatrixXd C_stud;
		kron(A, B, C_sol);
		test_kron(A, B, C_stud);
		CHECK(C_sol.norm() == doctest::Approx(C_stud.norm()).epsilon(1e-6));
	}
	
	SUBCASE("kron_mult(A,B)*x") {
		Eigen::VectorXd y_sol;
		Eigen::VectorXd y_stud;
		kron_mult(A, B, x, y_sol);
		test_kron_mult(A, B, x, y_stud);
		CHECK(y_sol.norm() == doctest::Approx(y_stud.norm()).epsilon(1e-6));
	}
	
	SUBCASE("kron_reshape(A,B)*x") {
		Eigen::VectorXd y_sol;
		Eigen::VectorXd y_stud;
		kron_reshape(A, B, x, y_sol);
		test_kron_reshape(A, B, x, y_stud);
		CHECK(y_sol.norm() == doctest::Approx(y_stud.norm()).epsilon(1e-6));
	}
}
