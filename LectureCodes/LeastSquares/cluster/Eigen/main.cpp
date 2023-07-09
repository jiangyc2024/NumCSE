///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "cluster.hpp"
#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXi;

int main () {
	
	// test cases
	Eigen::MatrixXd X1(2,4);
	X1 	<<	0,	4,	0,	4,
			0,	0,	1,	1;
	const int n1 = 2;
	Eigen::MatrixXd X2(2,10);
	X2 	<<	0,	4,	0,	4, 3, 5, 1, 0, 3, 7,
			0,	0,	1,	1, 4, 8, 3, 9, 2, 7;
	const int n2 = 3;
	Eigen::MatrixXd X3(3,6);
	X3 	<<	0,	4,	1,	4,	0,	1,
			0,	0,	0,	1,	4,	3,
			0,	0,	0,	0,	2,	2;
	const int n3 = 3;
	Eigen::MatrixXd X4(3,10);
	X4 	<<	0,	4,	0,	4, 3, 5, 1, 0, 3, 7,
			0,	0,	1,	1, 4, 8, 3, 9, 2, 7,
			3,	5,	3,	0, 4, 5, 9, 3, 5, 5;
	const int n4 =3;
	
	std::vector<Eigen::MatrixXd> X_vec = {X1, X2, X3, X4};
	std::vector<int> n_vec = {n1, n2, n3, n4};
	
	// solutions from former matlab code
	Eigen::MatrixXd C1(2,2);
	C1 	<<	0,    4.0000,
			0.5000,    0.5000;
	Eigen::VectorXi idx1(4);
	idx1 <<	1,     2,     1,     2;
	
	Eigen::MatrixXd C2(2,3);
	C2 	<<	3.5, 4, 1./3,
			1.75, 8, 4/3.;
	Eigen::VectorXi idx2(10);
	idx2 <<	 3,     1,     3,     1,     1,     2,     3,     2,     1,     2;
	
	Eigen::MatrixXd C3(3,3);
	C3 	<<	4,0.5,0.5,
			0.5,3.5,0,
			0,2,0;
	Eigen::VectorXi idx3(6);
	idx3 <<	 3,     1,     3,     1,     2,     2;
	
	Eigen::MatrixXd C4(3,3);
	C4 	<<	4./3,4,2.75,
			2/3.,8,2.25,
			2,13/3.,5.75;
	Eigen::VectorXi idx4(10);
	idx4 <<	1,     3,     1,     1,     3,     2,     3,     2,     3,     2;	
	
	std::vector<Eigen::MatrixXd> C_sol_vec = {C1, C2, C3, C4};
	std::vector<Eigen::VectorXi> idx_sol_vec = {idx1, idx2, idx3, idx4};

	for(unsigned int i = 0; i < X_vec.size(); ++i){
		std::cout << "##########  Test no. " << i + 1 << "  ##########" << std::endl;
		const Eigen::MatrixXd X = X_vec[i];
		const int n = n_vec[i];
		MatrixXd C;
		VectorXi idx;
		std::tie(C,idx) = cluster::pointcluster(X,n);
		std::cout << "C++, C matrix:\n" << C << "\nMatlab: C matrix:\n" << C_sol_vec[i] << std::endl;
		std::cout << "C++, idx vector:\n" << idx.transpose() << "\nMatlab: idx vector(C++ indices):\n" << (idx_sol_vec[i] - VectorXi::Ones(idx_sol_vec[i].size())).transpose() << std::endl << std::endl;
	}
	
	
	
	return 0;
}
