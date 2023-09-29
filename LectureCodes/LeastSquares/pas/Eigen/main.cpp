///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>

#include "pas.hpp"


int main () {
	Eigen::MatrixXd X(2,5);
	X 	<<	10,	2,	3,	4,	5,
			2,	4,	8,	1,	25;
	Eigen::VectorXi i1;
	Eigen::VectorXi i2;
	std::tie(i1, i2) = princaxissep::princaxissep(X);
	std::cout << i1 << std::endl << std::endl << i2 << std::endl;
	return 0;
}
