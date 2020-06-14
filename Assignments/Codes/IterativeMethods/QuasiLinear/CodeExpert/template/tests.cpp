#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

struct TestData {
	TestData() {
		n = 5;
		
		x = Eigen::VectorXd::Constant(n,1/std::sqrt(n));
		b = Eigen::VectorXd::Ones(n);
		
		atol = 1e-13;
		rtol = 1e-11;
	}
	int n;
	
    Eigen::VectorXd x;
    Eigen::VectorXd b;
    
    double atol;
    double rtol;
};

TestData data;

TEST_SUITE("QuasiLinear") {
	TEST_CASE("Eigen::VectorXd fixed_point_step" * doctest::description("Fixed point step")) {
		Eigen::VectorXd sol = fixed_point_step(data.x, data.b);
		Eigen::VectorXd stud = fixed_point_step_TEST(data.x, data.b);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("Eigen::VectorXd newton_step" * doctest::description("Newton step")) {
		Eigen::VectorXd sol = newton_step(data.x, data.b);
		Eigen::VectorXd stud = newton_step_TEST(data.x, data.b);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("Eigen::VectorXd solveQuasiLinSystem" * doctest::description("Fixed point iteration")) {
		Eigen::VectorXd sol = solveQuasiLinSystem(data.rtol, data.atol, data.b);
		Eigen::VectorXd stud = solveQuasiLinSystem_TEST(data.rtol, data.atol, data.b);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("Eigen::VectorXd solveQLSystem_Newton" * doctest::description("Newton iteration")) {
		Eigen::VectorXd sol = solveQLSystem_Newton(data.rtol, data.atol, data.b);
		Eigen::VectorXd stud = solveQLSystem_Newton_TEST(data.rtol, data.atol, data.b);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
}

