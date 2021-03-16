#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() { X = {0.0122, 0.3211, 0.7252, 1.3924, 1.5832}; }

  std::vector<double> X;
};

TestData data;

TEST_SUITE("CodeQuiz") {
  TEST_CASE("double myfunction_modified" *
            doctest::description("Test equivalent function")) {
    for (int i = 0; i < data.X.size(); i++) {
      double sol = myfunction_modified(data.X[i]);
      double stud = myfunction_modified_TEST(data.X[i]);

      CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
    }
  }

  TEST_CASE("double myfunction" * doctest::description("Unknown function") *
            doctest::skip()) {}
}
