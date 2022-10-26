#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <cmath>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("ChebPolyProperties") {
  TEST_CASE("bool checkDiscreteOrthogonality" *
            doctest::description("Check orthogonality of polynomials")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("Eigen::VectorXd bestpolchebnodes" *
            doctest::description("Interpolation")) {
    auto f = [](double x) { return 1. / (std::pow(5. * x, 2) + 1.); };
    Eigen::VectorXd sol = bestpolchebnodes(f, 12, 12);
    Eigen::VectorXd stud = bestpolchebnodes_TEST(f, 12, 12);
    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));

    sol = bestpolchebnodes(f, 20, 20);
    stud = bestpolchebnodes_TEST(f, 20, 20);
    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));

    sol = bestpolchebnodes(f, 20, 18);
    stud = bestpolchebnodes_TEST(f, 20, 18);
    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("std::vector<double> testBestPolyChebNodes" *
            doctest::description("Testing function")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
