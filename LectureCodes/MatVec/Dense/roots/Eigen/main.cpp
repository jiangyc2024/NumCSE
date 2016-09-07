#include <iostream>

#include <Eigen/Dense>

#include "zerosquadpol.hpp"

int main () {
	double alpha = 1, beta = -1;
	std::cout << zerosquadpol(alpha, beta) << std::endl;
	return 0;
}
