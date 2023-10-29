///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>


#include "cholfac.hpp"

int main () {
	const int n = 3;
	Eigen::MatrixXd A(n,n);
	Eigen::MatrixXd R(n,n);
	A << 	4, 12, -16,
			12, 37, -43,
			-16, -43, 98;
	cholfac::cholfac(A,R);
	std::cout << "A=\n" << A << std::endl;
	std::cout << "R=\n" << R << std::endl;
	std::cout << "R'*R=\n" << R.transpose()*R << std::endl;
	return 0;
}
