#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>
#include <vector>

struct TestData {
	TestData() {
		Psi = [] (double h, const Eigen::VectorXd & y0) -> Eigen::VectorXd { 
			return y0 * (1 + h); 
		};
		
		y0.resize(1);
		y0 << 1.0;
		T = 1.0;
		N = 8;
	}
	
	std::function<Eigen::VectorXd (double, Eigen::VectorXd)> Psi;
	Eigen::VectorXd y0;
	double T;
	int N;
};

TestData data;
TEST_SUITE("ODEsolve") {
	TEST_CASE("Vector psitilde" * doctest::description("Test psitilde")) {
		Eigen::VectorXd sol = psitilde(data.Psi, 1, 0.1, data.y0);
		Eigen::VectorXd stud = psitilde_TEST(data.Psi, 1, 0.1, data.y0);
	
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("std::vector<Vector> odeintequi" * doctest::description("Test equidistant integration")) {
		std::vector<Eigen::VectorXd> sol_vec = odeintequi(data.Psi, data.T, data.y0, data.N);
		std::vector<Eigen::VectorXd> stud_vec = odeintequi_TEST(data.Psi, data.T, data.y0, data.N);
		
		bool correct_num_steps = stud_vec.size() == sol_vec.size();
		CHECK(correct_num_steps);
		
		if (correct_num_steps) {
			for (int i = 0; i < sol_vec.size(); ++i) {
				Eigen::VectorXd sol = sol_vec[i];
				Eigen::VectorXd stud = stud_vec[i];
				
				bool vecsize_correct = sol.size() == stud.size();
				CHECK(vecsize_correct);
				
				if (vecsize_correct) {
					CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
				}
      		}
		}
	}
	
	TEST_CASE("std::pair< std::vector<double>, std::vector<Vector> > odeintssctrl" * doctest::description("Test adaptive integration")) {
		std::vector<Eigen::VectorXd> sol_vec = odeintssctrl(data.Psi, data.T, data.y0, 0.01, 1, 10e-5, 10e-5, 10e-5).second;
		std::vector<Eigen::VectorXd> stud_vec = odeintssctrl_TEST(data.Psi, data.T, data.y0, 0.01, 1, 10e-5, 10e-5, 10e-5).second;
		
		bool correct_num_steps = stud_vec.size() == sol_vec.size();
		CHECK(correct_num_steps);
		
		if (correct_num_steps) {
			for (int i = 0; i < sol_vec.size(); ++i) {
				Eigen::VectorXd sol = sol_vec[i];
				Eigen::VectorXd stud = stud_vec[i];
			
				bool vecsize_correct = sol.size() == stud.size();
				CHECK(vecsize_correct);
					
				if (vecsize_correct) {
					CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
				}
      		}
		}
	}
	
	TEST_CASE("double testcvpExtrapolatedEuler" * doctest::description("testcvpExtrapolatedEuler")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
	
	TEST_CASE("void solveTangentIVP" * doctest::description("solveTangentIVP")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
}

