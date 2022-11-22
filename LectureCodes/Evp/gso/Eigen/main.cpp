#include "gso.hpp"
#include <Eigen/Dense>
#include <Eigen/QR>
#include <iostream>


int main()
{
	const Eigen::Index n = 6;

	const Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
	const Eigen::MatrixXd Q = gso::gso(A);

	// check with householderQr
	auto qr = A.householderQr();
	const Eigen::MatrixXd Q_ex = qr.householderQ();
	std::cout << "Error = " << (Q_ex - Q).norm() << std::endl;

	std::cout << Q << std::endl << std::endl;
	std::cout << Q_ex << std::endl;

}
