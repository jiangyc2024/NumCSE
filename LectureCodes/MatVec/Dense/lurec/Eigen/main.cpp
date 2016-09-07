#include <iostream>

#include <Eigen/Dense>

#include "lurec.hpp"
#include "lurecdriver.hpp"

using namespace std;
using namespace Eigen;

int main () {
	int n = 3;
	MatrixXd A(n,n);
	A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	VectorXd b = VectorXd::LinSpaced(n,0,2);
	std::cout << "A=\n" << A << std::endl;
	MatrixXd L(n,n);
	MatrixXd U(n,n);
	lurecdriver(A, L, U);
	std::cout << "LU-Decomposition of A: \nL=\n" << L << "\nU=\n" << U << std::endl;
	// check / solve system
	VectorXd z = L.lu().solve(b);
	VectorXd xown = U.lu().solve(z);
	VectorXd xex = A.lu().solve(b);
	std::cout << "Own result:\n" << xown << std::endl;
	std::cout << "Eigen result:\n" << xex << std::endl;
	
	return 0;
}
