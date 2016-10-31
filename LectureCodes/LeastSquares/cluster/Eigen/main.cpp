///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>

#include "cluster.hpp"


int main () {
	//~ Eigen::MatrixXd X(2,4);
	//~ X 	<<	0,	4,	0,	4,
			//~ 0,	0,	1,	1;
	Eigen::MatrixXd X(2,10);
	X 	<<	0,	4,	0,	4, 3, 5, 1, 0, 3, 7,
			0,	0,	1,	1, 4, 8, 3, 9, 2, 7;
	int n = 3;
	MatrixXd C;
	VectorXi idx;
	std::tie(C,idx) = pointcluster(X,n);
	std::cout << C << std::endl << idx << std::endl;
	
	return 0;
}
