#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

struct TestData {
  TestData() {
    x.resize(8);
    x << 0.7, 3.3, 5.6, 7.5, 6.4, 4.4, 0.3, -1.1;
    y.resize(8);
    y << 4.0, 4.7, 4.0, 1.3, -1.1, -3.0, -2.5, 1.3;
  }

  Eigen::VectorXd x, y;
};

TestData data;

TEST_SUITE("Circle Approximation") {
  TEST_CASE("inline Eigen::MatrixXd lsqSVD" *
            doctest::description("SVD least-squares") * doctest::skip()) {}

  TEST_CASE("inline Eigen::MatrixXd lsqHHR" *
            doctest::description("QR least-squares") * doctest::skip()) {}

  TEST_CASE("inline Eigen::MatrixXd lsqNRM" *
            doctest::description("normal equations least-squares") *
            doctest::skip()) {}

  TEST_CASE("Eigen::Vector3d circl_alg_fit" *
            doctest::description("Algebraic fit")) {
    const Eigen::Vector3d sol = circl_alg_fit(data.x, data.y);
    const Eigen::Vector3d stud = circl_alg_fit_TEST(data.x, data.y);

    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-10));
  }

  TEST_CASE("void circl_geo_fit_GN" *
            doctest::description("Geometric fit Gauss-Newton")) {
    Eigen::Vector3d sol = circl_alg_fit(data.x, data.y);
    Eigen::Vector3d stud = sol;
    circl_geo_fit_GN(data.x, data.y, sol);
    circl_geo_fit_GN_TEST(data.x, data.y, stud);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-10));
  }

  TEST_CASE("void circl_geo_fit_N" *
            doctest::description("Geometric fit Newton")) {
    Eigen::Vector3d sol = circl_alg_fit(data.x, data.y);
    Eigen::Vector3d stud = sol;
    circl_geo_fit_N(data.x, data.y, sol);
    circl_geo_fit_N_TEST(data.x, data.y, stud);
    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-10));
  }

  TEST_CASE("void compare_convergence" *
            doctest::description("convergence plot")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("Eigen::Vector3d circl_svd_fit" *
            doctest::description("Constrained fit SVD")) {
    const Eigen::Vector3d sol = circl_svd_fit(data.x, data.y);
    const Eigen::Vector3d stud = circl_svd_fit_TEST(data.x, data.y);

    CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-10));
  }

  TEST_CASE("void plot" * doctest::description("Solution plot")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
