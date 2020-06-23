#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
	TestData() {
		n = 50;
		
		l_A = Eigen::VectorXd::LinSpaced(n - 1, 1, n - 1);
		d_A = Eigen::VectorXd::LinSpaced(n, n, 2 * n - 1);
		u_A = Eigen::VectorXd::LinSpaced(n - 1, 2 * n, 3 * n - 2);
		l_B = Eigen::VectorXd::Random(n - 1);
		d_B = Eigen::VectorXd::Random(n);
		u_B = Eigen::VectorXd::Random(n - 1);
	}
	
	std::size_t n;
	Eigen::VectorXd d_A, l_A, u_A, d_B, l_B, u_B;
};

TestData data;

TEST_SUITE("TriDiagQR") {
	TEST_CASE("inline double sign" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("template <typename T> int sgn" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("inline Eigen::Matrix2d Givens" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("std::tuple<double, double, double> compGivensRotation" * doctest::description("Givens rotation")) {
		Eigen::MatrixXd a = Eigen::MatrixXd::Random(2, 10);
		std::tuple<double, double, double> params_sol, params_stud;
		for (std::size_t k = 0; k < a.cols(); ++k) {
			params_sol = compGivensRotation(a.col(k));
			params_stud = compGivensRotation_TEST(a.col(k));
			CHECK(std::get<0>(params_sol) == doctest::Approx(std::get<0>(params_stud)).epsilon(1e-6));
			CHECK(std::get<1>(params_sol) == doctest::Approx(std::get<1>(params_stud)).epsilon(1e-6));
			CHECK(std::get<2>(params_sol) == doctest::Approx(std::get<2>(params_stud)).epsilon(1e-6));
		}
	}
	
	TEST_CASE("TriDiagonalQR [OUT OF CLASS]" * doctest::description("Constructor")) {
		TriDiagonalMatrix A(data.d_A, data.l_A, data.u_A);
		TriDiagonalQR sol_A(A);
		TriDiagonalQR_TEST stud_A(A);
		Eigen::MatrixXd Q_sol, Q_stud, R_sol, R_stud;

		sol_A.getDense(Q_sol, R_sol);
		stud_A.getDense(Q_stud, R_stud);
		CHECK((Q_sol - Q_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK((R_sol - R_stud).norm() == doctest::Approx(0.).epsilon(1e-6));

		TriDiagonalMatrix B(data.d_B, data.l_B, data.u_B);
		TriDiagonalQR sol_B(B);
		TriDiagonalQR_TEST stud_B(B);
		sol_B.getDense(Q_sol, R_sol);
		stud_B.getDense(Q_stud, R_stud);
		CHECK((Q_sol - Q_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		CHECK((R_sol - R_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("Eigen::VectorXd applyQT [OUT OF CLASS]" * doctest::description("applyQT()")) {
		TriDiagonalMatrix A(data.d_A, data.l_A, data.u_A);
		TriDiagonalQR sol_A(A);
		TriDiagonalQR_TEST stud_A(A);
		Eigen::VectorXd x = Eigen::VectorXd::Random(data.n);
		Eigen::VectorXd sol = sol_A.applyQT(x);
		Eigen::VectorXd stud = stud_A.applyQT(x);
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));

		TriDiagonalMatrix B(data.d_B, data.l_B, data.u_B);
		TriDiagonalQR sol_B(B);
		TriDiagonalQR_TEST stud_B(B);
		sol = sol_B.applyQT(x);
		stud = stud_B.applyQT(x);
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("Eigen::VectorXd solve [OUT OF CLASS]" * doctest::description("solve()")) {
		TriDiagonalMatrix A(data.d_A, data.l_A, data.u_A);
		TriDiagonalQR sol_A(A);
		TriDiagonalQR_TEST stud_A(A);

		TriDiagonalMatrix B(data.d_B, data.l_B, data.u_B);
		TriDiagonalQR sol_B(B);
		TriDiagonalQR_TEST stud_B(B);

		SUBCASE("normal test") {
			Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(data.n, 1, data.n);
			Eigen::VectorXd sol = sol_A.solve(b);
			Eigen::VectorXd stud = stud_A.solve(b);
			CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));

			b = Eigen::VectorXd::Random(data.n);
			sol = sol_B.solve(b);
			stud = stud_B.solve(b);
			CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
		}
		SUBCASE("exception test") {
			TriDiagonalMatrix B2 = B;
			Eigen::VectorXd b = Eigen::VectorXd::Random(data.n);
			B2.d(0) = 0.;
			B2.u(0) = 0.;
			TriDiagonalQR_TEST stud(B2);
			bool throws_runtime_error = false;
			Eigen::VectorXd x;
			try {
				x = stud.solve(b);
			} catch (std::runtime_error& err) {
				throws_runtime_error = true;
			}
			CHECK(throws_runtime_error);
		}

	}
	
	TEST_CASE("void getDense [OUT OF CLASS]" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("unsigned int invit" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("unsigned int invit<TriDiagonalMatrix>" * doctest::description("invit()")) {
		TriDiagonalMatrix A(data.d_A, data.l_A, data.u_A);
		TriDiagonalQR sol_A(A);
		TriDiagonalQR_TEST stud_A(A);
		Eigen::VectorXd x_sol = Eigen::VectorXd::Random(data.n);
		x_sol *= 10;
		Eigen::VectorXd x_stud = x_sol;
		unsigned int N_sol = invit<TriDiagonalMatrix>(A, x_sol);
		unsigned int N_stud = invit<TriDiagonalMatrix>(A, x_stud);
		CHECK(N_sol == N_stud);
		x_sol = Eigen::VectorXd::Random(data.n);
		x_sol *= 10;
		x_stud = x_sol;
		N_sol = invit<TriDiagonalMatrix>(A, x_sol, 1E-10, 100);
		N_stud = invit<TriDiagonalMatrix>(A, x_stud, 1E-10, 100);
		CHECK(N_sol == N_stud);
	}
}
