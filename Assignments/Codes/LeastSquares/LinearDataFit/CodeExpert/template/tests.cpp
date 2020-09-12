#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

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
	
	TEST_CASE("MatrixXd make_A" * doctest::description("make_A")) {
		Eigen::MatrixXd A_sol = make_A(data.t);
		Eigen::MatrixXd A_stud = make_A_TEST(data.t);
		
		REQUIRE(A_sol.rows() == A_stud.rows());
		REQUIRE(A_sol.cols() == A_stud.cols());
		CHECK((A_sol - A_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("VectorXd data_fit_normal" * doctest::description("normal equation")) {
		Eigen::VectorXd gamma_stud = data_fit_normal_TEST(data.t, data.f);
		
		REQUIRE(data.gamma1.size() == gamma_stud.size());
		CHECK((data.gamma1 - gamma_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("VectorXd data_fit_qr" * doctest::description("qr")) {
		Eigen::VectorXd gamma_stud = data_fit_qr_TEST(data.t, data.f);
		
		REQUIRE(data.gamma2.size() == gamma_stud.size());
		CHECK((data.gamma2 - gamma_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("VectorXd fit_plot" * doctest::description("fit_plot")) {
		Eigen::VectorXd yl_sol = fit_plot(data.gamma1, data.tl);
		Eigen::VectorXd yl_stud = fit_plot_TEST(data.gamma1, data.tl);
		
		REQUIRE(yl_sol.size() == yl_stud.size());
		CHECK((yl_sol - yl_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		
		yl_sol = fit_plot(data.gamma2, data.tl);
		yl_stud = fit_plot_TEST(data.gamma2, data.tl);
		
		REQUIRE(yl_sol.size() == yl_stud.size());
		CHECK((yl_sol - yl_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("VectorXd error_plot" * doctest::description("error_plot")) {
		Eigen::VectorXd err_sol = error_plot(data.gamma1);
		Eigen::VectorXd err_stud = error_plot_TEST(data.gamma1);
		
		REQUIRE(err_sol.size() == err_stud.size());
		CHECK((err_sol - err_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		
		err_sol = error_plot(data.gamma2);
		err_stud = error_plot_TEST(data.gamma2);
		
		REQUIRE(err_sol.size() == err_stud.size());
		CHECK((err_sol - err_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}

}

