#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"
#include <Eigen/Dense>

TEST_CASE("kron(A,B)*x") {
	Eigen::MatrixXd A(2, 2);
	A << 1, 2, 3, 4;
	Eigen::MatrixXd B(2, 2);
	B << 5, 6, 7, 8;
	Eigen::MatrixXd C_sol;
	Eigen::MatrixXd C_stud;
	srand(5);
	Eigen::VectorXd x = Eigen::VectorXd::Random(4);
	Eigen::VectorXd y_sol;
	Eigen::VectorXd y_stud;
	kron(A, B, C_sol);
	test_kron(A, B, C_stud);
	CHECK(C_sol.norm() - C_stud.norm() <= 1e-6);
}
