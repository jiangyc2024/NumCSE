#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include "timer.h"
using namespace std;
using namespace Eigen;
#include "linsolvetest.hpp"
int main () {
	std::cout << gausstiming() << std::endl;
	return 0;
}
