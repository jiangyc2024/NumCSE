///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>

#include "expeval.hpp"

int main () {
	int n = 20;
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(21,-20, n);
	Eigen::MatrixXd res(n+1,4);
	for(int i = 0; i <= n; ++i){
		res(i,0) = x(i);
		res(i,1) = expeval::expeval(x(i));
		res(i,2) = std::exp(x(i));
		res(i,3) = std::abs(res(i,2)-res(i,1))/res(i,2);
		// printing
		std::cout << std::fixed << std::setprecision(0) << std::setw(5) << res(i,0) << std::setprecision(10)<<
		std::setw(25) << std::scientific << res(i,1) << std::setw(25) << res(i,2) << 
		std::setw(25) << std::setprecision(15)<< std::fixed << res(i,3) << std::endl;
	}
	return 0;
}
