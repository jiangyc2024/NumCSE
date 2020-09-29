#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
	TestData() {
		n = 8;
		Xstop = MatrixXd(2,n);
		Xstop << 1, 3, 3, 1, -1, -3, -3, -1,  
                        -3, -1, 1, 3, 3, 1, -1, -3;
                Xpriority = MatrixXd(2,n);
		Xpriority <<  0, 3, 0, -3, 0, 2.5, 0, -2.5,
                             -3, 0, 3, 0, -2.5, 0, 2.5, 0;

                P = MatrixXd(2,n);
                P << 0.23657, 1.35369, -0.13624, -1.33702,
                     0.0989619, 0.993235, -0.0735973, -1.11657,
                    -2.76114, -2.60103, 2.90403, 2.66831,
                    -2.44302, -2.04656, 2.31922, 2.20296;
                
	}
	
	unsigned int n;
	Eigen::MatrixXd Xstop, Xpriority, P;
};

TestData data;

TEST_SUITE("ShapeIdent") {
	
	TEST_CASE("MatrixXd shape_ident_matrix" * doctest::description("shape_ident_matrix()")) {

		Eigen::MatrixXd sol, stud;
		
		sol = shape_ident_matrix(data.Xstop);
		stud = shape_ident_matrix_TEST(data.Xstop);
		
		REQUIRE(sol.size() == stud.size());
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	TEST_CASE("double solve_lsq" * doctest::description("solve_lsq()")) {

                MatrixXd A_stud, A_sol;
		double sol, stud;
		
		sol = solve_lsq(data.Xstop, data.P, A_sol);
		stud = solve_lsq_TEST(data.Xstop, data.P, A_stud);
		
		REQUIRE(A_sol.size() == A_stud.size());
		CHECK((A_sol - A_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
                CHECK(sol - stud == doctest::Approx(0.).epsilon(1e-6));
	}
        TEST_CASE("enum Shape" * doctest::description("skipped") * doctest::skip()) {}
        TEST_CASE("Shape identify" * doctest::description("identify()")) {

                MatrixXd A_stud, A_sol;
                Shape s_stud, s_sol;

                s_sol = identify(data.Xstop, data.Xpriority, data.P, A_sol);
                s_stud = identify_TEST(data.Xstop, data.Xpriority, data.P, A_stud);

                REQUIRE(A_sol.size() == A_stud.size());
		CHECK((A_sol - A_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
                CHECK(s_sol == s_stud);
        }
	
}

