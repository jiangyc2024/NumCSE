#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>
#include <vector>

struct TestData {
	TestData() {
		K = {2, 3, 5, 10, 15, 20};
	}
	
	std::vector<int> K;
};

TestData data;

TEST_SUITE("ExtrapolationToPi") {
	TEST_CASE("double extrapolate_to_pi" * doctest::description("Extrapolate to k")) {
		for (int i = 0; i < data.K.size(); i++) {
			double sol = extrapolate_to_pi(data.K[i]);
			double stud = extrapolate_to_pi_TEST(data.K[i]);
			
			CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
		}
	}
	
	TEST_CASE("void plotExtrapolationError" * doctest::description("Plot error") * doctest::skip()) {}
}

