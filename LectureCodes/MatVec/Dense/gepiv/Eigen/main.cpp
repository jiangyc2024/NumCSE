#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <limits>
using namespace std;
using namespace Eigen;
#include "gepiv.hpp"

int main () {
	int n = 3;
	MatrixXd A(n,n);
	A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	VectorXd b = VectorXd::LinSpaced(n,0,2);
	std::cout << "A=\n" << A << std::endl << "b=\n"<< b << std::endl;
	VectorXd x(n);
	gepiv(A, b, x);
	std::cout << std::endl << "x=\n" << x << std::endl;
	return 0;
}
