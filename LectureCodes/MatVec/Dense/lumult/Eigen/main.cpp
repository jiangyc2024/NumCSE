#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>
using namespace std;
using namespace Eigen;
#include "lumult.hpp"

int main () {
	int n = 3;
	MatrixXd A(n,n);
	MatrixXd L(n,n);
	MatrixXd U(n,n);
	L << 	1, 0, 0,
			-1, 1, 0,
			-2, 0.25, 1;
	U << 	1, 0, 2,
			0, 4, 3,
			0, 0, 5.25;
	lumult(L,U,A);
	
	std::cout << "L=\n" << L << std::endl << "U=\n" << U << std::endl;
	std::cout << "A=\n" << A << std::endl;
	
	/* should be equal
	 A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	 */ 
	return 0;
}
