#include <iostream>
#include <Eigen/Dense>

#include "cholfac.hpp"

int main () {
	int n = 3;
	Eigen::MatrixXd A(n,n);
	Eigen::MatrixXd R(n,n);
	A << 	4, 12, -16,
			12, 37, -43,
			-16, -43, 98;
	cholfac(A,R);
	std::cout << "A=\n" << A << std::endl;
	std::cout << "R=\n" << R << std::endl;
	std::cout << "R'*R=\n" << R.transpose()*R << std::endl;
	return 0;
}
