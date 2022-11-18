#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "copy.hpp"
#include "doctest.h"

// includes for test data
#include <Eigen/Dense>

struct TestData {
  TestData() {
    t = Eigen::VectorXd::LinSpaced(10, 0.1, 1.);
    f.resize(10);
    f << 100., 34., 17., 12., 9., 6., 5., 4., 4., 2.;
    tl = Eigen::VectorXd::LinSpaced(91, 0.1, 1.);
    gamma1 = data_fit_normal(t, f);
    gamma2 = data_fit_qr(t, f);
  }

  Eigen::VectorXd t, f, tl, gamma1, gamma2;
};

TestData data;

TEST_SUITE("LinearDataFit") {
  TEST_CASE("Eigen::MatrixXd make_A" * doctest::description("make_A")) {
    Eigen::MatrixXd A_sol = make_A(data.t);
    Eigen::MatrixXd A_stud = make_A_TEST(data.t);

    REQUIRE(A_sol.rows() == A_stud.rows());
    REQUIRE(A_sol.cols() == A_stud.cols());
    CHECK((A_sol - A_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("Eigen::VectorXd data_fit_normal" *
            doctest::description("normal equation")) {
    Eigen::VectorXd gamma_stud = data_fit_normal_TEST(data.t, data.f);

    REQUIRE(data.gamma1.size() == gamma_stud.size());
    CHECK((data.gamma1 - gamma_stud).norm() ==
          doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("Eigen::VectorXd data_fit_qr" * doctest::description("qr")) {
    Eigen::VectorXd gamma_stud = data_fit_qr_TEST(data.t, data.f);

    REQUIRE(data.gamma2.size() == gamma_stud.size());
    CHECK((data.gamma2 - gamma_stud).norm() ==
          doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("Eigen::VectorXd fitted_function" *
            doctest::description("fitted_function")) {
    Eigen::VectorXd yl_sol = fitted_function(data.gamma1, data.tl);
    Eigen::VectorXd yl_stud = fitted_function_TEST(data.gamma1, data.tl);

    REQUIRE(yl_sol.size() == yl_stud.size());
    CHECK((yl_sol - yl_stud).norm() == doctest::Approx(0.).epsilon(1e-8));

    yl_sol = fitted_function(data.gamma2, data.tl);
    yl_stud = fitted_function_TEST(data.gamma2, data.tl);

    REQUIRE(yl_sol.size() == yl_stud.size());
    CHECK((yl_sol - yl_stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }

  TEST_CASE("Eigen::VectorXd fitting_error" *
            doctest::description("fitting_error")) {
    Eigen::VectorXd err_sol = fitting_error(data.gamma1);
    Eigen::VectorXd err_stud = fitting_error_TEST(data.gamma1);

    REQUIRE(err_sol.size() == err_stud.size());
    CHECK((err_sol - err_stud).norm() == doctest::Approx(0.).epsilon(1e-8));

    err_sol = fitting_error(data.gamma2);
    err_stud = fitting_error_TEST(data.gamma2);

    REQUIRE(err_sol.size() == err_stud.size());
    CHECK((err_sol - err_stud).norm() == doctest::Approx(0.).epsilon(1e-8));
  }
}
