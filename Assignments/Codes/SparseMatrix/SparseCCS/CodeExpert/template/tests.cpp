#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
	TestData() {
		A = Eigen::MatrixXd::Zero(6, 6);
		A <<  4, -1,  0, -1,  0,  0,
			-1,  4, -1,  0, -1,  0,
			0, -1,  4,  0,  0, -1,
			-1,  0,  0,  4, -1,  0,
			0, -1,  0, -1,  4, -1,
			0,  0, -1,  0, -1,  4;
	}
	
	Eigen::MatrixXd A;
};

TestData data;

TEST_SUITE("SparseCCS") {
	
	TEST_CASE("void CCS" * doctest::description("tests val, row_ind, col_ptr")) {
		Eigen::VectorXd val_stud, val_sol, row_ind_stud, row_ind_sol, col_ptr_sol, col_ptr_stud;
		
		CCS(data.A, val_sol, row_ind_sol, col_ptr_sol);
		CCS_TEST(data.A, val_stud, row_ind_stud, col_ptr_stud);
		
		CHECK((val_sol - val_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK((row_ind_sol - row_ind_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK((col_ptr_sol - col_ptr_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}

}

