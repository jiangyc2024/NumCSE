#include <iostream>
#include <Eigen/Dense>
#include "levinson.hpp"

int main()
{
	Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(4, 0, 1);
	Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(10, 0, 1);
	Eigen::VectorXd x, y;

	levinson(u, b, x, y);
	std::cout << "x: " << std::endl << x << std::endl;
	std::cout << "y: " << std::endl << y << std::endl;
	return 0;
}
