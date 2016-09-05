#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include "timer.h"
using namespace std;
using namespace Eigen;
#include "linsolvetest.hpp"
int main () {
	std::cout << timing() << std::endl;
	int n = 10;
	MatrixXd A = MatrixXd::Random(n,n);
	VectorXd b = VectorXd::Random(n);
	VectorXd x(A.cols());
/* SAM_LISTING_BEGIN_0 */
// A is lower triangular
x = A.triangularView<Lower>().solve(b);
// A is upper triangular
x = A.triangularView<Upper>().solve(b);
// A is hermitian / self adjoint and positive definite
x = A.selfadjointView<Upper>().llt().solve(b);
// A is hermiatin / self adjoint (poitive or negative semidefinite)
x = A.selfadjointView<Upper>().ldlt().solve(b);
/* SAM_LISTING_END_0 */
	
	return 0;
}
