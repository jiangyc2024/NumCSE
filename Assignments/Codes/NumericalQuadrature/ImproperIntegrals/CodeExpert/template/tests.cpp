#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

struct TestData {
	TestData() {
		n = 10;
		
		f = [] (double t) {
			return 1 / (1 + std::pow(t,2));
		};
	}
	
	int n;
	std::function<double (double)> f;
};

TestData data;

TEST_SUITE("ImproperIntegrals") {
	TEST_CASE("double quadinf" * doctest::description("quadinf")) {
		double sol = quadinf(data.n, data.f);
		double stud = quadinf_TEST(data.n, data.f);
		
		CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("double quad" * doctest::description("Optional helper function")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}

	TEST_CASE("void cvgQuadInf" * doctest::description("Convergence")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
}

