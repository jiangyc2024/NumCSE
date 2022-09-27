#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <complex>

struct TestData {
  std::complex<double> w;
};

TestData data;

TEST_SUITE("ComplexRoot") {
  TEST_CASE("std::complex<double> myroot" *
            doctest::description("ComplexRoot")) {
    constexpr double eps = 1e-8;

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
