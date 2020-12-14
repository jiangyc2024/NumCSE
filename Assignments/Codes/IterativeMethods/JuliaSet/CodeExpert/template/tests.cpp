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
	TEST_CASE("Vector2d F" * doctest::description("Test F")) {
		Eigen::Vector2d sol = F(data.z);
		Eigen::Vector2d stud = F_TEST(data.z);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("Matrix2d DF" * doctest::description("Test DF")) {
		Eigen::Matrix2d sol = DF(data.z);
		Eigen::Matrix2d stud = DF_TEST(data.z);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("void julia" * doctest::description("Visualisation") * doctest::skip()) {}
}

