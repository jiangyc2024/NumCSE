#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

TEST_SUITE ("RecursionOrder") {
TEST_CASE("double testOrder" * doctest::description("testOrder")) {
      CHECK(testOrder(10) == doctest::Approx(testOrder_TEST(10)).epsilon(1e-9));
}
}

