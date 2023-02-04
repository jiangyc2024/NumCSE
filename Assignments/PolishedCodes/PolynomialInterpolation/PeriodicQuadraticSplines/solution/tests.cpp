#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <utility>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("PeriodicQuadraticSplines") {
  TEST_CASE("ClosedQuadraticSplineCurve [OUT OF CLASS]" *
            doctest::description("skipped") * doctest::skip()) {}

  TEST_CASE("bool checkC1 [OUT OF CLASS]" * doctest::description("skipped") *
            doctest::skip()) {}

  TEST_CASE("Eigen::Matrix2Xd curve_points [OUT OF CLASS]" *
            doctest::description("Testing output, not efficiency")) {
    // Comparing result with solution
    // Test case 1
    constexpr unsigned int n = 11;
    Eigen::MatrixXd p = Eigen::MatrixXd::Random(2, n);
    ClosedQuadraticSplineCurve_TEST obj(p);
    ClosedQuadraticSplineCurve objsol(p);

    Eigen::VectorXd v = Eigen::VectorXd::LinSpaced(15, 0, .99);
    Eigen::MatrixXd res = obj.curve_points(v);
    Eigen::MatrixXd ressol = objsol.curve_points(v);

    CHECK((res - ressol).norm() == doctest::Approx(0).epsilon(1e-10));

    // Test case 2
    constexpr unsigned int n1 = 21;
    Eigen::MatrixXd p1 = Eigen::MatrixXd::Random(2, n1);
    ClosedQuadraticSplineCurve_TEST obj1(p1);
    ClosedQuadraticSplineCurve objsol1(p1);

    Eigen::VectorXd v1 = Eigen::VectorXd::LinSpaced(25, 0, .99);
    Eigen::MatrixXd res1 = obj1.curve_points(v1);
    Eigen::MatrixXd ressol1 = objsol1.curve_points(v1);

    CHECK((res1 - ressol1).norm() == doctest::Approx(0).epsilon(1e-10));
  }

  TEST_CASE("Eigen::VectorXd local_curvatures [OUT OF CLASS]" *
            doctest::description("skipped") * doctest::skip()) {}

  TEST_CASE("double length [OUT OF CLASS]" * doctest::description("skipped") *
            doctest::skip()) {}
}
