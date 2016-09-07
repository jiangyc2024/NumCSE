#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include "apprpistable.hpp"

int main () {
	std::cout << setprecision(15) << std::fixed;
	Eigen::MatrixXd res = apprpistable();
	std::cout << res << std::endl;
	return 0;
}
