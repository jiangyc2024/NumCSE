#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include "ApproxPIinstable.hpp"

int main () {
	std::cout << setprecision(15) << std::fixed;
	Eigen::MatrixXd res = ApproxPIinstable();
	std::cout << res << std::endl;
	return 0;
}
