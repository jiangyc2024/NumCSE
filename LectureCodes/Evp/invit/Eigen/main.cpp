#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "invit.hpp"

int main()
{
	int n = 10;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
	A = A*A.transpose();

	double lmin;
	Eigen::VectorXd y;
	std::tie(lmin, y) = invit(A);

	std::cout << "lmin = " << lmin << std::endl;
	std::cout << "y = " << std::endl << y << std::endl;
	std::cout << " error norm(A*y - lmin*y) = " << (A*y - lmin*y).norm() << std::endl;

	// compare with eigensolver
	std::cout << "EigenSolver eigenvector = " << std::endl;
	Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
	int minIndex;
	solver.eigenvalues().real().minCoeff(&minIndex);
	std::cout << solver.eigenvectors().col(minIndex).real() << std::endl;
}
