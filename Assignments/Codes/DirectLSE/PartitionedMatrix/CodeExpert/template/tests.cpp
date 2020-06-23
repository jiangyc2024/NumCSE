#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
	TestData() {
		const std::size_t n = 50;
		R = Eigen::MatrixXd::Random(n, n).triangularView<Eigen::Upper>();
		v = Eigen::VectorXd::Random(n);
		u = Eigen::VectorXd::Random(n);
		bb = Eigen::VectorXd::Random(n + 1);
	}
	
	Eigen::MatrixXd R;
	Eigen::VectorXd v, u, bb;
};

TestData data;

TEST_SUITE("PartitionedMatrix") {
	
	TEST_CASE("void solvelse" * doctest::description("solvelse()")) {
		
		Eigen::VectorXd xo_sol, xo_stud;
		
		solvelse(data.R, data.v, data.u, data.bb, xo_sol);
		solvelse_TEST(data.R, data.v, data.u, data.bb, xo_stud);
		
		CHECK((xo_sol - xo_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		
	}
	
	TEST_CASE("bool testSolveLSE" * doctest::description("testSolveLSE")) {
		
		Eigen::VectorXd xe_sol, xe_stud;
		
		bool check_sol = testSolveLSE(data.R, data.v, data.u, data.bb, xe_sol);
		bool check_stud = testSolveLSE_TEST(data.R, data.v, data.u, data.bb, xe_stud);
		
		// check if LU-Decomposition was done correctly
		CHECK((xe_sol - xe_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK(check_stud == check_sol);
		
	}
	
}

