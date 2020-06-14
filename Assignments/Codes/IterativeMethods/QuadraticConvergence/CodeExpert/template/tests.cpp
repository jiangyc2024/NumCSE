#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

struct TestData {
	TestData() {
		x0 = 1.0;
		
		f = [] (double x) {
			return x * std::exp(x) - 1;
		}
	}
	
	double x0;
	std::function<double (double)> f;
};

TestData data;

TEST_SUITE("Steffensen") {
	TEST_CASE("template <class Function> double steffensen"* doctest::description("Testing Steffensen")) {
		// Test on the same function that testSteffensen() uses
		double sol = steffensen(data.f, data.x0);
		double stud = steffensen_TEST(data.f, data.x0);
		
		CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("void testSteffensen" * doctest::description("Test Steffensen") * doctest::skip()) {}
	
	TEST_CASE("template <class Function> double steffensen_log" * doctest::description("Steffensen log") * doctest::skip()) {}
	
	TEST_CASE("void orderSteffensen" * doctest::description("Steffensen order") * doctest::skip()) {}
}

