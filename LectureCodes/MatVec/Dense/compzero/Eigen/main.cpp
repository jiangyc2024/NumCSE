#include <iostream>
#include <cmath>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <figure/figure.hpp>
#include "../../roots/Eigen/zerosquadpol.hpp"
#include "../../zeroquadpolstab/Eigen/zerosquadpolstab.hpp"
#include "compzeros.hpp"

int main () {
	compzeros();
	return 0;
}
