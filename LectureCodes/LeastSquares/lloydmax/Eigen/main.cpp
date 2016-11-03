///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>

#include "lloydmax.hpp"


int main () {
	Eigen::MatrixXd X(2,5);
	X 	<<	10,	2,	3,	4,	5,
			2,	4,	8,	1,	25;
	Eigen::MatrixXd C(2,3);
	Eigen::VectorXi idx;
	Eigen::VectorXd cds;
	double tol = 1e-5;
	
	lloydmax(X, C, idx, cds, tol);
	std::cout << C << std::endl << idx << std::endl;
	
	return 0;
}
