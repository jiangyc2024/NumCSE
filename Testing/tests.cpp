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
	TestData() {
		// Initialise variables here
	}
	
	// Declare variables here
};

TestData data;

TEST_SUITE("_____EXERCISE_NAME_____") {
	// One test case for every function in the header files the student
	// is asked to modify.
	// 
	// FUNC of each TEST_CASE must be the part of the function's signature:
	// "RETURN_TYPE FUNCTION_NAME", e.g. "void kron". Don't include anything to do
	// with templates.
	//
	// Description can be anything. Usually taken from tests.csv or main.cpp.
	
	// Normal tests
	TEST_CASE("FUNCTION" * doctest::description("DESCRIPTION")) {
		// Student functions have _TEST appended to them.
		// Use CHECK macro to test. Semantics are CHECK(expr), pass if expr returns 1/true.
		// For floating point equality CHECK(x == doctest::Approx(0.).epsilon(1e-6));
		// The example here ensures compilation with a size check. Otherwise testing aborts
		// on run if the student's answer is unfinished or wrong.
		
		Eigen::VectorXd sol = FUNCTION();
		Eigen::VectorXd stud = FUNCTION_TEST();
		
		bool samesize = sol.size() == stud.size();
		CHECK(samesize);
		if (samesize) {
			CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		}
	}
	
	// Functions that print, but aren't tested.
	// Typically congerence/error studies, tabulating results.
	TEST_CASE("FUNCTION" * doctest::description("DESCRIPTION")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
	
	// Functions that aren't tested
	TEST_CASE("FUNCTION" * doctest::description("DESCRIPTION") * doctest::skip()) {}
}

