#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
};

TestData data;

TEST_SUITE("Cancellation") {
	
	TEST_CASE("void sinederv" * doctest::description("skipped") * doctest::skip()) {}
	
}

