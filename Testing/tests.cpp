//
// TEMPLATE FILE FOR INDEPENDENT TESTS OF HOMEWORK ASSIGNMENTS
// compatible with CodeExpert
// 
// Remove this after writing test
//

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
	// put all test data in here if all tests use the same data
	TestData() {
		// Initialise variables here
	}
	
	// Declare variables here
};

TestData data;

TEST_SUITE("EXERCISE_NAME" /* to be changed */) {
	
	// there must be a test case for every function in the header files that the
	// students should modify
	// the name of the TEST_CASE must match the convention "RETURN_TYPE FUNCTION_NAME"
	// e.g. "void kron", description can be anything
	// put a * doctest::skip() behind functions that should not be tested (as below)
	// If the function is templated, ignore the template line above the actual function
	// signature. 
	
	TEST_CASE("FUNCTION_SIG" /* to be changed */ * doctest::description("ANYTHING" /* to be changed */)) {
		
		// use doctest macros to test, e.g. CHECK(x <= y)
		// for floating point equality CHECK(x == doctest::Approx(0.).epsilon(1e-6))
		
	}
	
	// example for functions to be skipped by testing
	TEST_CASE("FUNCTION_SIG" /* to be changed */ * doctest::description("ANYTHING" /* to be changed */) * doctest::skip()) {}
}

