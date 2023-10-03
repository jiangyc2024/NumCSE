///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>

#include "arrowsystiming.hpp"

int main () {
	std::cout << std::scientific << std::setprecision(3) << arrowsys::arrowsystiming() << std::endl;
	return 0;
}
