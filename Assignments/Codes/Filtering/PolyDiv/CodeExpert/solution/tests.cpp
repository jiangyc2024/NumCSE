#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

struct TestData {
	TestData() {
        u.resize(4);
		v.resize(3);
	    u << 1, 2, 3, 4;
	    v << 10, 20, 30;
	}
    Eigen::VectorXd u;
	Eigen::VectorXd v;
};

TestData data;

TEST_SUITE("polyDiv") {
	
	TEST_CASE("VectorXd polyMult_naive" * doctest::description("Polynomial multiplication -- naive implementation")) {		
		auto sol = polyMult_naive(data.u, data.v);
		auto stud = polyMult_naive_TEST(data.u, data.v);
		
        CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}

	TEST_CASE("VectorXd polyMult_fast" * doctest::description("Polynomial multiplication -- efficient implementation")) {		
		auto sol = polyMult_fast(data.u, data.v);
		auto stud = polyMult_fast_TEST(data.u, data.v);
		
        CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
		
	TEST_CASE("VectorXd polyDiv" * doctest::description("Polynomial division -- efficient implementation")) {		
		auto sol_uv = polyMult_fast(data.u, data.v);
		auto stud_uv = polyMult_fast_TEST(data.u, data.v);
		auto sol_v = polyDiv(sol_uv, data.v);
		auto stud_v = polyDiv_TEST(stud_uv, data.v);
		
        CHECK((sol_v - stud_v).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK((data.v - stud_v).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
}
