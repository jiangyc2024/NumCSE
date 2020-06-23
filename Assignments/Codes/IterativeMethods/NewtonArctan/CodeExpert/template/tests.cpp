#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>
#include <vector>

struct TestData {
	TestData() {
		X = {0.21, 0.72, 1.3, 2.0, 3.5};
	}
	
	std::vector<double> X;
};

TestData data;

TEST_SUITE("NewtonArctan") {
	TEST_CASE("double newton_arctan" * doctest::description("arctan")) {
		for (int i = 0; i < data.X.size(); i++) {
			double sol = newton_arctan(data.X[i]);
			double stud = newton_arctan_TEST(data.X[i]);
			
			CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
		}
	}
}

