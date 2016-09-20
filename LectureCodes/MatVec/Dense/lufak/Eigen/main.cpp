///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <Eigen/Dense>

#include "lufak.hpp"

int main () {
	int n = 3;
	Eigen::MatrixXd A(n,n);
	Eigen::MatrixXd L(n,n);
	Eigen::MatrixXd U(n,n);
	A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	lufak(A,L,U);
	std::cout << "L=\n" << L << std::endl << "U=\n" << U << std::endl;
	return 0;
}
