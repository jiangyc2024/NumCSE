#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include "../../gausselimsolve/Eigen/gausselimsolve.hpp"
#include "../../lufak/Eigen/lufak.hpp"

int main () {
	/// tag
	#pragma begin<0>
	double epsilon = 5.0e-17;
	MatrixXd A(2,2), D(2,2);
	A << 	epsilon,	1.0,
			1.0,		1.0;
	D << 	1./epsilon,	0.0,
			0.0,		1.0;	
	VectorXd b(2), x2(2);
	b << 1.0, 2.0;
	A = D * A; 		b = D * b;
	VectorXd x1 = A.fullPivLu().solve(b);
	gausselimsolve(A,b, x2);
	MatrixXd L(2,2), U(2,2);
	lufak(A, L, U);
	VectorXd z = L.lu().solve(b);
	VectorXd x3 = U.lu().solve(z);
	cout << "x1=\n" << x1 << "\nx2\n" 
	<< x2 << "\nx3\n" << x3 << std::endl;
	#pragma end<0>
	return 0;
}