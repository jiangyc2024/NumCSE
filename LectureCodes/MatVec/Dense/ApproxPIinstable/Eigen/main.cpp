#include <iostream>
#include <Eigen/Dense>
#include <iomanip>
using namespace std;
using namespace Eigen;
#include "ApproxPIinstable.hpp"

int main () {
	std::cout << setprecision(15) << std::fixed;
	MatrixXd res = ApproxPIinstable();
	std::cout << res << std::endl;
	return 0;
}
