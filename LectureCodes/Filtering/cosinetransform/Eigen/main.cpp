#include "cosinetransform.hpp"

#include <Eigen/Dense>
#include <iostream>

int main()
{
	const Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(9, 0, 1);
	Eigen::VectorXd c;

	cosinetransform(y, c);
	std::cout << c << std::endl;
	return 0;
}
