#include <iostream>
#include <Eigen/Dense>
#include <cmath>
using namespace std;
using namespace Eigen;
#include "cholfac.hpp"

int main () {
	int n = 3;
	MatrixXd A(n,n);
	MatrixXd R(n,n);
	A << 	4, 12, -16,
			12, 37, -43,
			-16, -43, 98;
	cholfac(A,R);
	std::cout << "A=\n" << A << std::endl;
	std::cout << "R=\n" << R << std::endl;
	std::cout << "R'*R=\n" << R.transpose()*R << std::endl;
	return 0;
}
