///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <limits>

using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "gsrecpiv.hpp"
/*
 * wrong/error code! --> look at lupiv
 */

//NOLINTBEGIN(bugprone-exception-escape)
int main () {
	const int n = 3;
	MatrixXd A(n,n);
	A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	const VectorXd b = VectorXd::LinSpaced(n,0,2);
	std::cout << "A=\n" << A << std::endl;
	MatrixXd L(n,n);
	MatrixXd U(n,n);
	MatrixXd Adec = gsrecpiv::gsrecpiv(A);
	std::cout << "A decomposition=\n" << Adec << std::endl;
	U = Adec.triangularView<Eigen::Upper>();
	L.setIdentity();
	L += Adec.triangularView<Eigen::StrictlyLower>();
	std::cout << "LU-Decomposition of A: \nL=\n" << L << "\nU=\n" << U << std::endl;
	// check / solve system
	const VectorXd z = L.lu().solve(b);
	const VectorXd xown = U.lu().solve(z);
	const VectorXd xex = A.lu().solve(b);
	std::cout << "Own result:\n" << xown << std::endl;
	std::cout << "Eigen result:\n" << xex << std::endl;
	
	return 0;
}
//NOLINTEND(bugprone-exception-escape)