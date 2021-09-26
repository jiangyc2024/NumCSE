#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <complex>

struct TestData {
	std::complex<double> w;
};

TestData data;

TEST_SUITE("ComplexRoot") {
	
	TEST_CASE("std::complex<double> myroot" * doctest::description("ComplexRoot")) {
		const double eps = 1e-6;
		
		data.w = std::complex<double>(1e20, 5);
		auto sol = myroot(data.w);
		auto stud = myroot_TEST(data.w);
		CHECK(sol.real() == doctest::Approx(stud.real()).epsilon(eps));
		CHECK(sol.imag() == doctest::Approx(stud.imag()).epsilon(eps));
		
		data.w = std::complex<double>(-5, 1e20);
		sol = myroot(data.w);
		stud = myroot_TEST(data.w);
		CHECK(sol.real() == doctest::Approx(stud.real()).epsilon(eps));
		CHECK(sol.imag() == doctest::Approx(stud.imag()).epsilon(eps));
		
		data.w = std::complex<double>(1e-8, 0);
		sol = myroot(data.w);
		stud = myroot_TEST(data.w);
		CHECK(sol.real() == doctest::Approx(stud.real()).epsilon(eps));
		CHECK(sol.imag() == doctest::Approx(stud.imag()).epsilon(eps));
	}

}

