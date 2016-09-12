#include <iostream>

#include <Eigen/Dense>

#include "lupiv.hpp"
#include "lupivdriver.hpp"

int main () {
	int n = 3;
	Eigen::MatrixXd A(n,n);
	A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(n,0,2);
	std::cout << "A=\n" << A << std::endl;
	Eigen::MatrixXd L(n,n);
	Eigen::MatrixXd U(n,n);
	lupivdriver(A,L,U);
	std::cout << "LU-Decomposition of A: \nL=\n" << L << "\nU=\n" << U << std::endl;
	// check / solve system
	Eigen::VectorXd z = L.lu().solve(b);
	Eigen::VectorXd xown = U.lu().solve(z);
	Eigen::VectorXd xex = A.lu().solve(b);
	std::cout << "Own result:\n" << xown << std::endl;
	std::cout << "Eigen result:\n" << xex << std::endl;
	
	return 0;
}
