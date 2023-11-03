///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
///            R. Hiptmair <hiptmair@sam.math.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>


using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "smw.hpp"

int main () {
	const int n = 10;
	const MatrixXd A = MatrixXd::Random(n,n);
	const VectorXd b = VectorXd::Random(n);
	const VectorXd u = VectorXd::Random(n);
	const VectorXd v = VectorXd::Random(n);
	const Eigen::PartialPivLU<MatrixXd> lu(A);
	std::cout << "Solving rank-1 updated LSE" << std::endl;
	std::cout << smw::smw(u,v,lu,b) << std::endl << std::endl;
	std::cout << (A + u * v.transpose()).lu().solve(b) << std::endl;
	return 0;
}

/* OLD VERSION 
int main () {
	int n = 10;
	MatrixXd A = MatrixXd::Random(n,n);
	VectorXd b = VectorXd::Random(n);
	VectorXd u = VectorXd::Random(n);
	VectorXd v = VectorXd::Random(n);
	Eigen::PartialPivLU<MatrixXd> lu(A);
	MatrixXd LU = lu.matrixLU();
	MatrixXd L = MatrixXd::Identity(n,n);
	L.triangularView<StrictlyLower>() = LU;
	MatrixXd U = LU.triangularView<Upper>();
	MatrixXd P = lu.permutationP();
	// Permutation adjustment is bypassed
	std::cout << smw(L,U,u,v,b) << std::endl << std::endl;
	std::cout << (L*U + u * v.transpose()).lu().solve(b) << std::endl;
	return 0;
}
*/
