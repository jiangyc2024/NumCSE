#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
	TestData() {
		M = 100;
		xa = VectorXd::Random(M);
	}
	
	unsigned int M;
	Eigen::VectorXd xa;
};

TestData data;

TEST_SUITE("StructuredMatrixVector") {

	TEST_CASE("void multAminSlow" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("void multAmin" * doctest::description("multAmin")) {
		Eigen::VectorXd sol, stud;
		
		multAmin(data.xa, sol);
		multAmin_TEST(data.xa, stud);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("int multAmin_runtime" * doctest::description("multAmin_runtime")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
	
	TEST_CASE("Eigen::MatrixXd multABunitv" * doctest::description("multABunitv")) {
		Eigen::MatrixXd sol, stud;
		
		sol = multABunitv();
		stud = multABunitv_TEST();
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
}

