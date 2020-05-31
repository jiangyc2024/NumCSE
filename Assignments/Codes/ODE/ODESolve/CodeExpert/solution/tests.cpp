#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>
#include <vector>

struct TestData {
	TestData() {
		Psi = [] (double h, const Eigen::VectorXd & y0) -> Eigen::VectorXd { 
			return (1+h)*y0; 
		};
		
		y0 << 1.0;
		T = 1.0;
		N = 8;
	}
	
	auto Psi;
	Eigen::VectorXd y0(1);
	double T;
	int N;
};

TestData data;

TEST_SUITE("ODEsolve") {
	// TODO: these might not work because they're templated
	TEST_CASE("Vector psitilde" * doctest::description("Test psitilde") * doctest::skip()) {
		Eigen::VectorXd sol = psitilde(data.Psi, 1, 0.1, data.y0);
		Eigen::VectorXd stud = psitilde_TEST(data.Psi, 1, 0.1, data.y0);
	
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("std::vector<Vector> odeintequi" * doctest::description("Test equidistant integration")) {
		std::vector<Eigen::VectorXd> sol = odeintequi(data.Psi, data.T, data.y0, data.N);
		std::vector<Eigen::VectorXd> stud = odeintequi_TEST(data.Psi, data.T, data.y0, data.N);
		
		for (int i = 0; i < data.N; ++i) {
			CHECK((sol[i](0) - stud[i](0)).norm() == doctest::Approx(0.).epsilon(1e-6));
      	}
	}
	
	// TODO: this one additionally has a multiline signature
	TEST_CASE("std::pair< std::vector<double>, std::vector<Vector> > odeintssctrl" * doctest::description("Test adaptive integration")) {
		std::vector<Eigen::VectorXd> sol = odeintssctrl(data.Psi, data.T, data.y0, 0.01, 1, 10e-5, 10e-5, 10e-5).second;
		std::vector<Eigen::VectorXd> stud = odeintssctrl_TEST(data.Psi, data.T, data.y0, 0.01, 1, 10e-5, 10e-5, 10e-5).second;
	
		for (int i = 0; i < data.N; ++i) {
			CHECK((sol[i](0) - stud[i](0)).norm() == doctest::Approx(0.).epsilon(1e-6));
      	}
	}
	
	TEST_CASE("double testcvpExtrapolatedEuler" * doctest::description("testcvpExtrapolatedEuler") * doctest::skip()) {}
	
	TEST_CASE("void solveTangentIVP" * doctest::description("solveTangentIVP") * doctest::skip()) {}
}

