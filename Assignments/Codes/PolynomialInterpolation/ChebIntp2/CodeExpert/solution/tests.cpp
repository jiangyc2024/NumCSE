#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

struct TestData {
	TestData() {
		a = 0;
		b = 1;
		N = 100;
		tol = 1e-4;
		
		f = [] (double t) {
			return std::sin(std::exp(2 * t));
		};
	}
	
	double a;
	double b;
	unsigned int N;
	double tol;
	
	std::function<double (double)> f;
};

TestData data;
TEST_SUITE("AdaptivePolyIntp") {
	TEST_CASE("VectorXd adaptivepolyintp" * doctest::description("Nodes")) {
		Eigen::VectorXd sol_vec = adaptivepolyintp(data.f, data.a, data.b, data.tol, data.N);
		Eigen::VectorXd stud_vec = adaptivepolyintp_TEST(data.f, data.a, data.b, data.tol, data.N);
		
		bool samesize = sol_vec.size() == stud_vec.size();
		CHECK(samesize);
		
		if (samesize) {
			for (int i = 0; i < sol_vec.size(); i++) {
				double sol = sol_vec[i];
				double stud = stud_vec[i];
				
				CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
			}
		}
	}
	
	TEST_CASE("void plotInterpolationError" * doctest::description("Plot error") * doctest::skip()) {}
}

