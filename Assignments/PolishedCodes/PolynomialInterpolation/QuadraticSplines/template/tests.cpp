#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    t.resize(3);
    t << .1, .2, .5;

    y.resize(4);
    y << 1, 2, 1, 2;

    x << 0, 0.4, 0.9;
  }
  Eigen::VectorXd t;
  Eigen::VectorXd y;
  Eigen::Vector3d x;
};

TestData data;
TEST_SUITE("QuadraticSplines") {
  TEST_CASE("std::pair<Eigen::VectorXd, Eigen::VectorXd> increments" *
            doctest::description("increments")) {
    std::pair<Eigen::VectorXd, Eigen::VectorXd> sol = increments(data.t);
    std::pair<Eigen::VectorXd, Eigen::VectorXd> stud = increments_TEST(data.t);

    const bool samesize = sol.first.size() == stud.first.size() &&
                          sol.second.size() == stud.second.size();
    REQUIRE(samesize);
    CHECK((sol.first - stud.first).norm() == doctest::Approx(0.).epsilon(1e-8));
    CHECK((sol.second - stud.second).norm() ==
          doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("Eigen::VectorXd compute_c" * doctest::description("compute_c")) {
    Eigen::VectorXd sol = compute_c(data.t, data.y);
    Eigen::VectorXd stud = compute_c_TEST(data.t, data.y);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("Eigen::VectorXd compute_d" * doctest::description("compute_d")) {
    Eigen::VectorXd c = compute_c(data.t, data.y);

    Eigen::VectorXd sol = compute_d(c, data.t);
    Eigen::VectorXd stud = compute_d_TEST(c, data.t);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("Eigen::VectorXd quadspline" * doctest::description("fval")) {
    Eigen::VectorXd sol = quadspline(data.t, data.y, data.x);
    Eigen::VectorXd stud = quadspline_TEST(data.t, data.y, data.x);

    REQUIRE(sol.size() == stud.size());
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("void plotquadspline" *
            doctest::description("Plot quadtratic spline")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("std::vector<double> qsp_error" *
            doctest::description("Error of quadtratic spline")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
