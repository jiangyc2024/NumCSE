///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>

#include "lupiv.hpp"
#include "lupivdriver.hpp"

//NOLINTBEGIN(bugprone-exception-escape)
int main () {
	const int n = 3;
	Eigen::MatrixXd A(n,n);
	A << 	1, 0, 2,
			-1, 4, 1,
			-2, 1, 2;
	const Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(n,0,2);
	std::cout << "A=\n" << A << std::endl;
	Eigen::MatrixXd L(n,n);
	Eigen::MatrixXd U(n,n);
	lupiv::lupivdriver(A,L,U);
	std::cout << "LU-Decomposition of A: \nL=\n" << L << "\nU=\n" << U << std::endl;
	// check / solve system
	const Eigen::VectorXd z = L.lu().solve(b);
	const Eigen::VectorXd xown = U.lu().solve(z);
	const Eigen::VectorXd xex = A.lu().solve(b);
	std::cout << "Own result:\n" << xown << std::endl;
	std::cout << "Eigen result:\n" << xex << std::endl;
	
	return 0;
}
//NOLINTEND(bugprone-exception-escape)