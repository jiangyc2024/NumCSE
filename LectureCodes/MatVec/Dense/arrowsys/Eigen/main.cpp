#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include <vector>
#include "timer.h"
using namespace std;
using namespace Eigen;
#include "arrowsys_fast.hpp"
#include "arrowsys_slow.hpp"
#include "arrowsys_sparse.hpp"
#include "arrowsystiming.hpp"
int main () {
	std::cout << arrowsystiming() << std::endl;
	return 0;
}
