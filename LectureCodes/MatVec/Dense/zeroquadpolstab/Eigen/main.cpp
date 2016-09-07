#include <iostream>

#include <Eigen/Dense>

#include "zerosquadpolstab.hpp"

int main () {
	double alpha = 1, beta = -1;
	std::cout << zerosquadpolstab(alpha, beta) << std::endl;
	return 0;
}
