#include "eigqr.hpp"
#include <Eigen/Eigenvalues>
#include <iostream>


int main()
{
	const Eigen::Index n = 6;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
	A = A*A.transpose();
	auto ev = eigqr::eigqr(A);
	
	// compare with EigenSolver
	const Eigen::EigenSolver<Eigen::MatrixXd> eig(A, false);
	auto const & ev2 = eig.eigenvalues();

	std::cout << "QR-algorithm:" << std::endl << ev << std::endl << std::endl;
	std::cout << "EigenSolver" << std::endl << ev2 << std::endl << std::endl;
}
