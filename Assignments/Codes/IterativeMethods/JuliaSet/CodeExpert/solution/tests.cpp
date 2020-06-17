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
	TEST_CASE("Vector2d F_vector" * doctest::description("Test F")) {
		Eigen::Vector2d sol = F_vector(data.z);
		Eigen::Vector2d stud = F_vector_TEST(data.z);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("Matrix2d DF_matrix" * doctest::description("Test DF")) {
		Eigen::Matrix2d sol = DF_matrix(data.z);
		Eigen::Matrix2d stud = DF_matrix_TEST(data.z);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("void julia" * doctest::description("Visualisation") * doctest::skip()) {}
}

