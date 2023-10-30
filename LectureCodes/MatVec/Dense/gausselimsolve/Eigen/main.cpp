///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>

#include "gausselimsolve.hpp"

int main () {
	const int n = 3;
	Eigen::MatrixXd A(n,n);
	
	A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	const Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(n,0,2);
	std::cout << "A=\n" << A << std::endl << "b=\n"<< b << std::endl;
	Eigen::VectorXd x(n);
	gausselimsolve::gausselimsolve(A, b, x);
	std::cout << "x=\n" << x << std::endl;
	return 0;
}
