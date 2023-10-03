///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>


#include <Eigen/Dense>

#include "ApproxPIinstable.hpp"

int main () {
	std::cout << std::setprecision(15) << std::fixed;
	const Eigen::MatrixXd res = ApproxPIinstable::ApproxPIinstable();
	std::cout << res << std::endl;
	return 0;
}
