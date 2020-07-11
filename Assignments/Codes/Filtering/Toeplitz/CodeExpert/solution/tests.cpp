#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

TEST_SUITE("Toeplitz") {

	TEST_CASE("MatrixXd toeplitz" * doctest::description("Toeplitz matrix")) {
		Eigen::VectorXd c = Eigen::VectorXd::LinSpaced(50, 1, 50);
		Eigen::VectorXd r = Eigen::VectorXd::LinSpaced(50, 1, 40);

		auto sol = toeplitz(c, r);
		auto stud = toeplitz_TEST(c, r);

		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("VectorXd toepmatmult" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("VectorXd toepmult" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("VectorXd ttmatsolve" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("VectorXd ttrecsolve" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("VectorXd ttsolve" * doctest::description("ttsolve wrapper")) {
		Eigen::VectorXd h = Eigen::VectorXd::LinSpaced(50, 1, 50);
		Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(50, 51, 100);

		auto sol = ttsolve(h, y);
		auto stud = ttsolve_TEST(h, y);

		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
}

