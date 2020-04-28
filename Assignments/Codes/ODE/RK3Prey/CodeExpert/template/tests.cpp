#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"
#include <Eigen/Dense>

struct TestData {
};

TestData data;

TEST_SUITE("RK3Prey") {
	TEST_CASE("RKIntegrator" * doctest::description("...") * doctest::skip()) {}
	
	TEST_CASE("std::vector<State> solve" * doctest::description("...") * doctest::skip()) {}
	
	TEST_CASE("void step" * doctest::description("...") * doctest::skip()) {}
	
	TEST_CASE("double RK3prey" * doctest::description("...") * doctest::skip()) {}
}
