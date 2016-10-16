///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <numeric>

#include <Eigen/Dense>

#include "tridiagqr.hpp"

int main () {
	std::cout << "#####################################" << std::endl;
	std::cout << "Test of tridiagqr" << std::endl;
	std::cout << "#####################################" << std::endl;
	Eigen::VectorXd c(3),d(4),e(3),b(4);
	c << 1,2,5;
	d << 4,6,5,9;
	e << 5,2,4;
	/*
		[	4,	1,	0,	0;
			5,	6,	2,	0;
			0,	2,	5	5;
			0,	0,	4,	9]
 	 */
 	b << 1, 2, 3, 4;
	std::cout << tridiagqr(c,d,e,b) << std::endl;
	return 0;
}
