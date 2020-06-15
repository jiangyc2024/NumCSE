#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

using namespace Eigen;

struct TestData {
	TestData() {
		t.resize(3);
		t << .1, .2, .5;
		
		y.resize(4);
		y << 1, 2, 1, 2;
		
		x << 0, 0.4, 0.9;
	}
	VectorXd t;
	VectorXd y;
	Vector3d x;
};

TestData data;
TEST_SUITE("QuadraticSplines") {
	TEST_CASE("std::pair<VectorXd, VectorXd> increments" * doctest::description("increments")) {
		std::pair<VectorXd, VectorXd> sol = increments(data.t);
		std::pair<VectorXd, VectorXd> stud = increments_TEST(data.t);
		
		CHECK((sol.first - stud.first).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK((sol.second - stud.second).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("VectorXd compute_c" * doctest::description("compute_c")) {
		VectorXd sol = compute_c(data.t, data.y);
		VectorXd stud = compute_c_TEST(data.t, data.y);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("VectorXd compute_d" * doctest::description("compute_d")) {
		VectorXd c = compute_c(data.t, data.y);
		
		VectorXd sol = compute_d(c, data.t);
		VectorXd stud = compute_d_TEST(c, data.t);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("VectorXd quadspline" * doctest::description("fval")) {
		VectorXd sol = quadspline(data.t, data.y, data.x);
		VectorXd stud = quadspline_TEST(data.t, data.y, data.x);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("void plotquadspline" * doctest::description("Plot quadtratic spline") * doctest::skip()) {}
	
	TEST_CASE("std::vector<double> qsp_error" * doctest::description("Error of quadtratic spline") * doctest::skip()) {}
}

