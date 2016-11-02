#include <Eigen/Dense>
#include <Eigen/QR>
#include <iostream>
#include "gso.hpp"


int main()
{
	int n = 4;

	Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
	Eigen::MatrixXd Q = gso(A);

	// check if orthonormal
	std::cout << Q*Q.transpose() << std::endl;
}
