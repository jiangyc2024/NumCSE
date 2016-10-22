#include <iostream>
#include <Eigen/Dense>
#include "sspowitrpex.hpp"

int main()
{
	int n = 10;
	Eigen::VectorXd d = Eigen::VectorXd::LinSpaced(n,1,n);
	sspowitrpex(d);
}


