#include <iostream>
#include <Eigen/Dense>
#include "sspowitrp.hpp"


int main()
{
	// Generate matrix
	int n = 10;
	Eigen::VectorXd d = Eigen::VectorXd::LinSpaced(n,1,n);
	Eigen::VectorXd oneToN = Eigen::VectorXd::LinSpaced(n,1,n);
	Eigen::MatrixXd Z = Eigen::MatrixXd::Ones(n,n);
	Z += oneToN.array().sqrt().matrix().asDiagonal();

	auto qr = Z.householderQr();   // generate orthogonal matrix
	Eigen::MatrixXd Q = qr.householderQ();
	Eigen::MatrixXd A = Q*d.asDiagonal()*Q.transpose(); // ``synthetic'' \Blue{$\VA=\VA^T$} with spectrum \Blue{$\sigma(\VA) =\{d_1,\ldots,d_n\}$}

	Eigen::MatrixXd V;
	Eigen::VectorXd ev;

	// call power iteration
	int k = 4;
	sspowitrp(A, k, 8, 1e-6, 100, ev, V);

	std::cout << std::endl << std::endl;
	std::cout << "Largest 4 eigenvalues:" << std::endl;
	std::cout << ev << std::endl << std::endl;
	std::cout << "Largest 4 eigenvectors:" << std::endl;
	std::cout << V << std::endl;
}
