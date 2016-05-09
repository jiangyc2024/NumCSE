#include <iostream>
#include <Eigen/Dense>
#include <cmath>
using namespace std;
using namespace Eigen;
#include "blockgs.hpp"

int main () {
	int n = 3;
	MatrixXd A(n,n);
	A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	VectorXd b = VectorXd::LinSpaced(n,0,2);
	std::cout << "A=\n" << A << std::endl << "b=\n"<< b << std::endl;
	blockgs(A);
	std::cout << "After gaussian elimination: A=\n" << A << std::endl;
	return 0;
}
