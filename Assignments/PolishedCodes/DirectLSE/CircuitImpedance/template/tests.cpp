#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

constexpr double eps = 1e-12;

TEST_SUITE("Resistance to impedance map") {
  TEST_CASE("NodalPotentials [OUT OF CLASS]" *
            doctest::description("NodalPotentials constructor") *
            doctest::skip()) {}

  TEST_CASE("Eigen::VectorXd operator() [OUT OF CLASS]" *
            doctest::description("NodalPotentails operator()")) {
    // all resistors have same resistance
    NodalPotentials sol1(1, 1);
    NodalPotentials_TEST stud1(1, 1);
    // check for two different source voltages
    Eigen::VectorXd vec_sol1 = sol1(1);
    Eigen::VectorXd vec_stud1 = stud1(1);
    REQUIRE(vec_sol1.size() == vec_stud1.size());
    CHECK((vec_sol1 - vec_stud1).norm() == doctest::Approx(0.).epsilon(eps));
    vec_sol1 = sol1(5.1);
    vec_stud1 = stud1(5.1);
    CHECK((vec_sol1 - vec_stud1).norm() == doctest::Approx(0.).epsilon(eps));

    // check for different R and Rx
    NodalPotentials sol2(2, 3.5);
    NodalPotentials_TEST stud2(2, 3.5);
    // check for different source voltages
    Eigen::VectorXd vec_sol2 = sol2(1);
    Eigen::VectorXd vec_stud2 = stud2(1);
    REQUIRE(vec_sol2.size() == vec_stud2.size());
    CHECK((vec_sol2 - vec_stud2).norm() == doctest::Approx(0.).epsilon(eps));
    vec_sol2 = sol2(5.1);
    vec_stud2 = stud2(5.1);
    CHECK((vec_sol2 - vec_stud2).norm() == doctest::Approx(0.).epsilon(eps));
  }

  TEST_CASE("ImpedanceMap [OUT OF CLASS]" *
            doctest::description("ImpedanceMap constructor") *
            doctest::skip()) {}

  TEST_CASE("double operator() [OUT OF CLASS]" *
            doctest::description("ImpedanceMap operator()")) {
    ImpedanceMap sol1(1, 1);
    ImpedanceMap_TEST stud1(1, 1);
    for (std::size_t i = 1; i <= 1024; i *= 2) {
      CHECK(sol1(i) == doctest::Approx(stud1(i)).epsilon(eps));
    }
    ImpedanceMap sol2(1, 5.1);
    ImpedanceMap_TEST stud2(1, 5.1);
    for (std::size_t i = 1; i <= 1024; i *= 2) {
      CHECK(sol2(i) == doctest::Approx(stud2(i)).epsilon(eps));
    }
  }
}
