///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "linsolvetest.hpp"

int main () {
	std::cout << timing() << std::endl;
	int n = 10;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
	Eigen::VectorXd b = Eigen::VectorXd::Random(n);
	Eigen::VectorXd x(A.cols());
/* SAM_LISTING_BEGIN_0 */
// A is lower triangular
x = A.triangularView<Eigen::Lower>().solve(b);
// A is upper triangular
x = A.triangularView<Eigen::Upper>().solve(b);
// A is hermitian / self adjoint and positive definite
x = A.selfadjointView<Eigen::Upper>().llt().solve(b);
// A is hermiatin / self adjoint (poitive or negative semidefinite)
x = A.selfadjointView<Eigen::Upper>().ldlt().solve(b);
/* SAM_LISTING_END_0 */
	
	return 0;
}
