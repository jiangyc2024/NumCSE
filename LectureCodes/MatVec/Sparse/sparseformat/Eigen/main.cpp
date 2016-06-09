#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <vector>
#include <limits>
#include <figure/figure.hpp>
#include "timer.hpp"
#include "spdiags.hpp"
using namespace std;
using namespace Eigen;

int main () {
	/*
	// Test
	/// Matlab example 5A
	VectorXi diag_no(3);
	diag_no << -2,0,2;
	MatrixXd B(5,3);
	B << 	1,	6,	11,
			2,	7,	12,
			3,	8,	13,
			4,	9,	14,
			5,	10,	15;
	std::cout << B << std::endl;
	std::cout << spdiags(B, diag_no, 5,5) << std::endl;
	std::cout << spdiags(B, diag_no, 5,4) << std::endl;
	std::cout << spdiags(B, diag_no, 4,5) << std::endl;
	*/
	#pragma begin<0>
	RowVectorXd diag_el(3);
	diag_el << -1, 2, 5;
	auto sparseM = [&](double n){	// lambda function for
		VectorXi diag_no(3);		// generation of matrix
		diag_no << -n/2, 0, n/2;
		MatrixXd B = diag_el.replicate(n,1);
		return spdiags(B, diag_no, n,n);
	};
	mgl::Figure fig1;
	fig1.spy(MatrixXd(sparseM(16)));	fig1.setFontSize(5);
	fig1.title("Pattern for matrix {\\bf A} for n = 16");
	fig1.save("spdiagsmatspy_cpp");
	// Timing
	int nruns = 3;
	double t1, t2; int min_i = 1, max_i = 20;
	MatrixXd times(max_i-min_i+1,3); 
	Timer timer;	// timer class
	double tmp = 0; // prevent optimization
	for(int i = min_i; i <= max_i ; ++i){
		int n = std::pow(2,i), m = n/2;
		SparseMatrix<double> A = sparseM(n);
		RowVectorXd v1(n); VectorXd v2(n);
		double t1 = std::numeric_limits<double>::max(); double t2 = t1;
		for(int k = 0; k < nruns; ++k){
			timer.start();
			for(int j = 1; j <= 5; ++j)
				v1 += A.row(m);
			timer.stop();
			t1 = std::min(t1, timer.duration());
			timer.start();
			for(int j = 1; j <= 5; ++j)
				v2 += A.col(m);
			timer.stop();
			t2 = std::min(t2, timer.duration());
		}
		times(i-min_i,0) = n; times(i-min_i,1) = t1;
		times(i-min_i,2) = t2;
		tmp += v1*v2; // prevent optimization
	}
	// plotting ...
	#pragma end<0>
	
	
	std::cout << tmp << std::endl;
	//std::cout << times << std::endl;
	mgl::Figure fig2;
	fig2.setlog(true,true);
	fig2.setFontSize(5);
	fig2.plot(times.col(0), times.col(1), " +r-").label("row access");
	fig2.plot(times.col(0), times.col(2), " *b-").label("column access");
	fig2.plot(times.col(0), times(0,1)*times.col(0)/times(0,0), " k-")
												.label("O(n)");
	fig2.xlabel("size n of sparse quadratic matrix");
	fig2.ylabel("access time [s]");
	fig2.legend(0.05,0.95);
	fig2.save("sparseaccess_cpp");
	return 0;
}
