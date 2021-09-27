#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("Cancellation") {
  TEST_CASE("void sinederv" * doctest::description("skipped")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
