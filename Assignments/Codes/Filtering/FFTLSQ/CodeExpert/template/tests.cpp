#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

TEST_SUITE("FFTLSQ") {
	TEST_CASE("VectorXd gauss_fit" * doctest::description("gauss_fit()")) {
		Eigen::VectorXd d(10);
		d << 0.987214,
		1.03579,
		0.997689,
		0.917471,
		1.00474,
		0.92209,
		1.03517,
		1.08863,
		0.904992,
		0.956089;
		
		auto sol = gauss_fit(d, 3);
		auto stud = gauss_fit_TEST(d, 3);
		
		REQUIRE(sol.size() == stud.size());
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
}


