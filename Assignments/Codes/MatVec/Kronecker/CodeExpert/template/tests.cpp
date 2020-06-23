#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"
#include <Eigen/Dense>

struct TestData {
	TestData() {
		A = Eigen::MatrixXd::Random(10, 10);
		B = Eigen::MatrixXd::Random(10, 10);
		x = Eigen::VectorXd::Random(100);
	}
	Eigen::MatrixXd A;
	Eigen::MatrixXd B;
	Eigen::VectorXd x;
};

TestData data;

TEST_SUITE("Kronecker") {
	TEST_CASE("void kron" * doctest::description("kron(A,B)*x")) {
		Eigen::MatrixXd C_sol;
		Eigen::MatrixXd C_stud;
		kron(data.A, data.B, C_sol);
		kron_TEST(data.A, data.B, C_stud);
		CHECK((C_sol - C_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("void kron_mult" * doctest::description("kron_mult(A,B)*x")) {
		Eigen::VectorXd y_sol;
		Eigen::VectorXd y_stud;
		kron_mult(data.A, data.B, data.x, y_sol);
		kron_mult_TEST(data.A, data.B, data.x, y_stud);
		CHECK((y_sol - y_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("void kron_reshape" * doctest::description("kron_reshape(A,B)*x")) {
		Eigen::VectorXd y_sol;
		Eigen::VectorXd y_stud;
		kron_reshape(data.A, data.B, data.x, y_sol);
		kron_reshape_TEST(data.A, data.B, data.x, y_stud);
		CHECK((y_sol - y_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("void kron_runtime" * doctest::description("Test runtime")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
}
