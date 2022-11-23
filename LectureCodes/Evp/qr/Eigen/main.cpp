#include "eigqr.hpp"
#include <Eigen/Eigenvalues>
#include <iostream>


int main()
{
	int n = 6;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
	A = A*A.transpose();
	auto ev = eigqr(A);
	
	// compare with EigenSolver
	Eigen::EigenSolver<Eigen::MatrixXd> eig(A, false);
	auto ev2 = eig.eigenvalues();

	std::cout << "QR-algorithm:" << std::endl << ev << std::endl << std::endl;
	std::cout << "EigenSolver" << std::endl << ev2 << std::endl << std::endl;
}
