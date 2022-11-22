#include "fftbasedsolutionlocal.hpp"
#include <Eigen/Dense>
#include <iostream>

int main()
{
	const Eigen::VectorXd b1 = Eigen::VectorXd::LinSpaced(4, 1, 4);
	const Eigen::VectorXd b2 = Eigen::VectorXd::LinSpaced(10, 1, 10);
	const Eigen::MatrixXd B = b1*b2.transpose();
	std::cout << B << std::endl << std::endl;

	const double c = 1.0;
	const double cx = 2.0;
	const double cy = 3.0;

	Eigen::MatrixXd X;
	fftbasedsolutionlocal(B, c, cx, cy, X);
	std::cout << X << std::endl << std::endl;

	return 0;
}
