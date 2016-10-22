#include <Eigen/Dense>
#include <Eigen/QR>
#include <iostream>
#include "gso.hpp"


int main()
{
	int n = 6;

	Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
	Eigen::MatrixXd Q = gso(A);

	// check with householderQr
	auto qr = A.householderQr();
	Eigen::MatrixXd Q_ex = qr.householderQ();
	std::cout << "Error = " << (Q_ex - Q).norm() << std::endl;

	std::cout << Q << std::endl << std::endl;
	std::cout << Q_ex << std::endl;

}
