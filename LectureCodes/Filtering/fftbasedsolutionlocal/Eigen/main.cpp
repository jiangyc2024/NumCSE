#include <iostream>
#include <Eigen/Dense>
#include "fftbasedsolutionlocal.hpp"

int main()
{
	Eigen::VectorXd b1 = Eigen::VectorXd::LinSpaced(4, 1, 4);
	Eigen::VectorXd b2 = Eigen::VectorXd::LinSpaced(10, 1, 10);
	Eigen::MatrixXd B = b1*b2.transpose();
	std::cout << B << std::endl << std::endl;

	double c, cx, cy;
	c = 1.0;
	cx = 2.0;
	cy = 3.0;

	Eigen::MatrixXd X;
	fftbasedsolutionlocal(B, c, cx, cy, X);
	std::cout << X << std::endl << std::endl;

	return 0;
}
