#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

using namespace Eigen;

struct TestData {
	TestData() {
		T = 1;
		n = 3;
		Y0.resize(n, n);
		Y0 << 1, 1, 0, 0, 3, 2, 1, 5, 2;
	}
	
	double T;
	unsigned int n;
	MatrixXd Y0;
};

TestData data;

TEST_SUITE("NLMatODE") {
	TEST_CASE("MatrixXd matode" * doctest::description("Check matode matrix at T")) {
		MatrixXd sol = matode(data.Y0, data.T);
		MatrixXd stud = matode_TEST(data.Y0, data.T);

		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}

	TEST_CASE("bool checkinvariant" * doctest::description("Test whether invariant for Y0 was preserved or not")) {
		bool sol = checkinvariant(data.Y0, data.T);
		bool stud = checkinvariant_TEST(data.Y0, data.T);
		
		CHECK(sol == stud);
	}
	
	TEST_CASE("double cvgDiscreteGradientMethod" * doctest::description("Test fitted convergance rate")) {
		double sol = cvgDiscreteGradientMethod();
		double stud = cvgDiscreteGradientMethod_TEST();
		
		CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
	}
}

