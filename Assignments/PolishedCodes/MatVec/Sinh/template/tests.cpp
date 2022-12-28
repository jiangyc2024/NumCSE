#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

TEST_SUITE("Sinh") {
  TEST_CASE("double sinh_unstable" * doctest::description("unstable sinh") *
            doctest::skip()) {}

  TEST_CASE("void sinhError" * doctest::description("Error tabulation")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
