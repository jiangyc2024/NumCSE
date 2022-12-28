#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    // use the image shot at focus 2
    double focus = 2;
    M = set_focus(focus);
  }
  Eigen::MatrixXd M;
};

TestData data;

TEST_SUITE("Autofocus") {
  TEST_CASE("Eigen::MatrixXd set_focus" * doctest::description("skipped") *
            doctest::skip()) {}

  TEST_CASE("double high_frequency_content" *
            doctest::description("Compute V")) {
    double sol = high_frequency_content(data.M);
    double stud = high_frequency_content_TEST(data.M);

    CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("double autofocus" *
            doctest::description("Find most focused image")) {
    double stud = autofocus_TEST();

    // opitmal focus is at around 2
    CHECK(std::abs(stud) == doctest::Approx(2.).epsilon(1e-1));
  }

  TEST_CASE("void save_image" * doctest::description("Save image")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("void plot_freq" * doctest::description("Plot spectrum")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("void plotV" * doctest::description("Plot V")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
