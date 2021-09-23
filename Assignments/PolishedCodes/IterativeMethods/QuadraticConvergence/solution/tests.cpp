#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

TEST_SUITE("Steffensen") {
	TEST_CASE("double steffensen" * doctest::description("Testing Steffensen")) {
    constexpr double x0 = 1.;
    auto f = [] (double x) { return x * std::exp(x) - 1; };
		// Test on the same paramters that testSteffensen() uses
		const double sol = steffensen(f, x0);
		const double stud = steffensen_TEST(f, x0);
		
		CHECK(sol == doctest::Approx(stud).epsilon(1e-9));
	}
	
	TEST_CASE("void testSteffensen" * doctest::description("Test Steffensen")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
	
	TEST_CASE("double steffensen_log" * doctest::description("Steffensen log")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
	
	TEST_CASE("void orderSteffensen" * doctest::description("Steffensen order")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
}

