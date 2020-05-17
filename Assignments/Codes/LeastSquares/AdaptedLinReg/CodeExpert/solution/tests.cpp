#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
	TestData() {
		n = 25;
		t = Eigen::VectorXd::LinSpaced(n, 0., 1.);
		std::srand(41);
		noise = Eigen::VectorXd::Random(n);
	}
	
	int n;
	Eigen::VectorXd t, noise;
};

TestData data;

TEST_SUITE("AdaptedLinReg") {
	
	TEST_CASE("VectorXd linReg" * doctest::description("linReg()")) {
		Eigen::VectorXd y = 12 * data.t + Eigen::VectorXd::Constant(data.n, -154) + data.noise;
		Eigen::VectorXd sol, stud;
		
		sol = linReg(data.t, y);
		stud = linReg_TEST(data.t, y);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("VectorXd expFit" * doctest::description("expFit()")) {
		Eigen::VectorXd y = 17 * Eigen::exp(-3 * data.t.array()).matrix() + 0.5 * data.noise;
		Eigen::VectorXd sol, stud;
		
		sol = expFit(data.t, y);
		stud = expFit_TEST(data.t, y);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
}

