#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

TEST_SUITE("MatrixClass") {
	
	TEST_CASE("Eigen::Matrix<double, 2, 2> smallTriangular" * doctest::description("smallTriangular")) {
		
		const auto sol = smallTriangular(1, 2, 3);
		const auto stud = smallTriangular_TEST(1, 2, 3);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		
	}
	
	TEST_CASE("Eigen::MatrixXd constantTriangular" * doctest::description("constantTriangular")) {
		
		const auto sol = constantTriangular(3, 20);
		const auto stud = constantTriangular_TEST(3, 20);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		
	}
	
	TEST_CASE("double casting" * doctest::description("casting")) {
		
		const auto sol = casting();
		const auto stud = casting_TEST();
		
		CHECK(sol == doctest::Approx(stud).epsilon(1e-6));
		
	}
	
	TEST_CASE("Eigen::VectorXcd arithmetics" * doctest::description("arithmetics")) {
		
		const auto sol = arithmetics(4);
		const auto stud = arithmetics_TEST(4);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		
	}

}

