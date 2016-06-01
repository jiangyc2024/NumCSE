#include <iostream>
#include <cmath>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include "zerosquadpolstab.hpp"

int main () {
	double alpha = 1, beta = -1;
	std::cout << zerosquadpolstab(alpha, beta) << std::endl;
	return 0;
}
