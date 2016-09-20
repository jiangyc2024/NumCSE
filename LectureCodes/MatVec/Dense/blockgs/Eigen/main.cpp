///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>

#include "blockgs.hpp"

int main () {
	int n = 3;
	Eigen::MatrixXd A(n,n);
	A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(n,0,2);
	Eigen::MatrixXd Ab(n,n+1);
	Ab << A, b;
	std::cout << "[A,b]=\n" << Ab << std::endl;
	blockgs(Ab);
	std::cout << "After gaussian elimination: Ab=\n" << Ab << std::endl;
	return 0;
}
