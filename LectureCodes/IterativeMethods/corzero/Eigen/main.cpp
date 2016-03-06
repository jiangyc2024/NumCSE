#include <Eigen/Dense>
#include <iostream>
#include "corzero.hpp"
#include "figure.hpp"

int main()
{
	Eigen::VectorXd rates;
	Eigen::VectorXd err;
	
	corzero(0.4, rates, err);

	return 0;
}
