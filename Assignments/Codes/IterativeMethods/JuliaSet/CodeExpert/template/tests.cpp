#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

struct TestData {
	TestData() {
		z.resize(2);
		z << -0.888, 1.333;
	}
	
	Eigen::Vector2d z;
};

TestData data;

TEST_SUITE("JuliaSet") {
	TEST_CASE("Vector2d juliaF" * doctest::description("Test F")) {
		Eigen::Vector2d sol = juliaF(data.z);
		Eigen::Vector2d stud = juliaF_TEST(data.z);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));
	}
	
	TEST_CASE("Matrix2d juliaDF" * doctest::description("Test DF")) {
		Eigen::Matrix2d sol = juliaDF(data.z);
		Eigen::Matrix2d stud = juliaDF_TEST(data.z);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-9));
	}
	
	TEST_CASE("void julia" * doctest::description("Visualisation")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}

