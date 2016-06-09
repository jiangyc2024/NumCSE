#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include "timer.hpp"
using namespace std;
using namespace Eigen;
#include "arrowsys_fast.hpp"
#include "arrowsys_slow.hpp"
#include "arrowsystiming.hpp"
int main () {
	std::cout << arrowsystiming() << std::endl;
	return 0;
}
