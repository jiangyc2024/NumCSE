#include "cosinetransform.hpp"
#include "icosinetransform.hpp"
#include <Eigen/Dense>
#include <iostream>

int main()
{
	Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(9, 0, 1);
	Eigen::VectorXd c;

	cosinetransform(y, c);
	icosinetransform(c, y);
	std::cout << y << std::endl;
	return 0;
}
