#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>
#include <vector>

struct TestData {
	TestData() {
		// Definitiion of coefficients in Butcher scheme
		s = 2;
		A.resize(s, s);
		b.resize(s);
	    A << 5. / 12., -1. / 12., 3. / 4., 1. / 4.;
        b << 3. / 4., 1. / 4.;
        
        // Coefficients and handle for prey/predator model
        double alpha1 = 3.;
        double alpha2 = 2.;
        double beta1 = 0.1;
        double beta2 = 0.1;
        
        f = [&alpha1, &alpha2, &beta1, &beta2](const VectorXd &y) {
			auto temp = y;
			temp(0) *= alpha1 - beta1 * y(1);
			temp(1) *= -alpha2 + beta2 * y(0);
			return temp;
		};
		
		Jf = [&alpha1, &alpha2, &beta1, &beta2](const VectorXd &y) {
			MatrixXd temp(2, 2);
			temp << alpha1 - beta1 * y(1), -beta1 * y(0), beta2 * y(1),
				-alpha2 + beta2 * y(0);
			return temp;
		};
		
		// Dimension of state space
 		d = 2;
 		
 		
 		// Initial value for model
 		y0.resize(d);
 		y0 << 100, 5;
 		
 		// Variables for testing (iterations, final times)
 		N = {128, 256, 512, 1024, 2048};
        T = {2.0, 5.0, 10.0, 20.0, 50.0};  
	}
	
	unsigned int s;
	MatrixXd A;
	VectorXd b;
	
	std::function<VectorXd(const VectorXd&)> f;
	std::function<MatrixXd(const VectorXd&)> Jf;
	
	unsigned int d;
	VectorXd y0;
	
	std::vector<unsigned int> N;
	std::vector<double> T;
};

TestData data;

TEST_SUITE("ImplRK Prey") {
	TEST_CASE("std::vector<VectorXd> solve" * doctest::description("RK solving")) {
		// Instantiate classes
		implicit_RKIntegrator RK_sol(data.A, data.b);
		implicit_RKIntegrator_TEST RK_stud(data.A, data.b);
		
		// Test with various values
		for (int i = 0; i < data.T.size(); i++) {
			for (int j = 0; j < data.N.size(); j++) {
				std::vector<VectorXd> sol_vec = RK_sol.solve(data.f, data.Jf, data.T[i], data.y0, data.N[j]);
				std::vector<VectorXd> stud_vec = RK_stud.solve_TEST(data.f, data.Jf, data.T[i], data.y0, data.N[j]);				
				
				VectorXd sol = sol_vec.back();
				VectorXd stud = stud_vec.back();
				
				
				CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
				
				
				/*
				// DEBUGGING
				// TODO: copy this tests.cpp to the solution one too
				// TODO: figure out why the tests trivially pass.
				// It looks like it's only looking at the first elem of the 
				// eigen vector
				// It's returning exactly the value of y0 each time
				// This is the default behaviour of the template file
				// It looks like the solution is identical to the
				// template.
				
				printf("\n\nSOL VECTOR OF TEST #%d: \n\n", (i + j));

				printf("{");
				for (int k = 0; k < sol_vec.size(); k++) {
					
					printf("(");
					for (int l = 0; l < sol_vec[k].size(); l++) {
						double elem = (sol_vec[k])[l];
						printf("%d, ", elem);
					}
					printf(")\n");
				}
				printf("\n}\n\n");
				*/
			}
		}
	}
	
	TEST_CASE("MatrixXd kron" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("implicit_RKIntegrator" * doctest::description("skipped") * doctest::skip()) {}
}

