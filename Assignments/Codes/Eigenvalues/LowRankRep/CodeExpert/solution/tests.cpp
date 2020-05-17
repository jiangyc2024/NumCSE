#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
	TestData() {
		m = 3;
		n = 2;
		k = 2;
		X.resize(m, n);
		X << 5, 0, 2, 1, 7, 4;
		A.resize(m, k);
		B.resize(n, k);
		A << 2, 1, 2, 3, 6, 1;
		B << 4, 4, 5, 0;
		Ax.resize(m, k);
		Ay.resize(m, k);
		Bx.resize(n, k);
		By.resize(n, k);
		Ax << 1, 0, 9, 2, 6, 3;
		Ay << 8, -2, 3, 4, 5, 8;
		Bx << 2, 1, 2, 3;
		By << 4, 4, 5, 0;
	}
	
	std::size_t m, n, k;
	Eigen::MatrixXd X, A, B, Ax, Ay, Bx, By;
};

TestData data;

TEST_SUITE("LowRankRep") {
	
	TEST_CASE("std::pair<MatrixXd,MatrixXd> factorize_X_AB" * doctest::description("factor")) {
		auto sol = factorize_X_AB(data.X, data.k);
		auto stud = factorize_X_AB_TEST(data.X, data.k);
		
		CHECK((sol.first - stud.first).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK((sol.second - stud.second).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("std::tuple<MatrixXd, MatrixXd, MatrixXd> svd_AB" * doctest::description("svd U, S, V")) {
		auto sol = svd_AB(data.A, data.B);
		auto stud = svd_AB_TEST(data.A, data.B);
		
		CHECK((std::get<0>(sol) - std::get<0>(stud)).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK((std::get<1>(sol) - std::get<1>(stud)).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK((std::get<2>(sol) - std::get<2>(stud)).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("std::pair<MatrixXd,MatrixXd> rank_k_approx" * doctest::description("approx Z")) {
		auto sol = rank_k_approx(data.Ax, data.Ay, data.Bx, data.By);
		auto stud = rank_k_approx_TEST(data.Ax, data.Ay, data.Bx, data.By);
		
		CHECK((sol.first - stud.first).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK((sol.second - stud.second).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
}

