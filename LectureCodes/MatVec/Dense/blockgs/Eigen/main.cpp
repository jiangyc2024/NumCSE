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
	MatrixXd Ab(n,n+1);
	Ab << A, b;
	std::cout << "[A,b]=\n" << Ab << std::endl;
	blockgs(Ab);
	std::cout << "After gaussian elimination: Ab=\n" << Ab << std::endl;
	return 0;
}
