#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip>
using namespace std;
using namespace Eigen;
#include "apprpistable.hpp"

int main () {
	std::cout << setprecision(15) << std::fixed;
	MatrixXd res = apprpistable();
	std::cout << res << std::endl;
	return 0;
}
