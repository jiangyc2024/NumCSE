#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "solution.hpp"
#include "periodicquadraticsplines.hpp"
#include <utility>
#include <cmath>


TEST_SUITE("CompactStorageQR") {

  TEST_CASE("curve_points()" * doctest::description("Testing output, not efficiency")) {
    // Comparing result with solution
    // Test case 1
    unsigned n = 11;
    Eigen::MatrixXd p = Eigen::MatrixXd::Random(2,n);
    ClosedQuadraticSplineCurve obj(p);
    solution::ClosedQuadraticSplineCurve objsol(p);

    Eigen::VectorXd v = Eigen::VectorXd::LinSpaced(15,0,.99);
    Eigen::MatrixXd res = obj.curve_points(v);
    Eigen::MatrixXd ressol = objsol.curve_points(v);

    CHECK((res-ressol).norm() == doctest::Approx(0).epsilon(1e-10));

    // Test case 2
    unsigned n1 = 21;
    Eigen::MatrixXd p1 = Eigen::MatrixXd::Random(2,n1);
    ClosedQuadraticSplineCurve obj1(p1);
    solution::ClosedQuadraticSplineCurve objsol1(p1);

    Eigen::VectorXd v1 = Eigen::VectorXd::LinSpaced(25,0,.99);
    Eigen::MatrixXd res1 = obj1.curve_points(v1);
    Eigen::MatrixXd ressol1 = objsol1.curve_points(v1);

    CHECK((res1-ressol1).norm() == doctest::Approx(0).epsilon(1e-10));
  }

}
