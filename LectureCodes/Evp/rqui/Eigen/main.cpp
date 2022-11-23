#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <iostream>
#include "rqui.hpp"


int main()
{
	int n = 10;
	Eigen::VectorXd d = Eigen::VectorXd::LinSpaced(n,1,n);

	// This is for producing matrix A from vector d
	Eigen::MatrixXd Z = d.array().sqrt().matrix().asDiagonal();
	Z += Eigen::MatrixXd::Ones(n,n);
	Eigen::MatrixXd Q = Z.colPivHouseholderQr().matrixQ();

	Eigen::MatrixXd A = Q*d.asDiagonal()*Q.transpose();
	A = A*A.transpose();
	Eigen::SparseMatrix<double> As = A.sparseView();

	double lmin;
	Eigen::VectorXd y;
	std::tie(y, lmin) = rqui(As);

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
