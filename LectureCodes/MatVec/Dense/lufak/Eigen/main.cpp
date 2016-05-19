#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>
using namespace std;
using namespace Eigen;
#include "lufak.hpp"

int main () {
	int n = 3;
	MatrixXd A(n,n);
	MatrixXd L(n,n);
	MatrixXd U(n,n);
	A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	lufak(A,L,U);
	std::cout << "L=\n" << L << std::endl << "U=\n" << U << std::endl;
	return 0;
}
