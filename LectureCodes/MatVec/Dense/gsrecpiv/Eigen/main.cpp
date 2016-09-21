///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <limits>
using namespace std;
using namespace Eigen;
#include "gsrecpiv.hpp"
/*
 * wrong/error code! --> look at lupiv
 */

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
	MatrixXd Adec = gsrecpiv(A);
	std::cout << "A decomposition=\n" << Adec << std::endl;
	U = Adec.triangularView<Upper>();
	L.setIdentity();
	L += Adec.triangularView<StrictlyLower>();
	std::cout << "LU-Decomposition of A: \nL=\n" << L << "\nU=\n" << U << std::endl;
	// check / solve system
	VectorXd z = L.lu().solve(b);
	VectorXd xown = U.lu().solve(z);
	VectorXd xex = A.lu().solve(b);
	std::cout << "Own result:\n" << xown << std::endl;
	std::cout << "Eigen result:\n" << xex << std::endl;
	
	return 0;
}
