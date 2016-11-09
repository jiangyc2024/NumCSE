#include <iostream>
#include <Eigen/Dense>
#include "cosinetransform.hpp"

int main()
{
	Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(9, 0, 1);
	Eigen::VectorXd c;

	cosinetransform(y, c);
	std::cout << c << std::endl;
	return 0;
}
