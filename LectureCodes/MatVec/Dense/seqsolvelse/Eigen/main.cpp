///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>

#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;

int main () {
	const int N = 10;
	const int n = 100;
	VectorXd x(n);
	x.setOnes();
	MatrixXd A(n,n);
	A.setRandom();
	auto some_function = [](VectorXd &x){return x*42;};
	VectorXd b = some_function(x);
	
	
/* SAM_LISTING_BEGIN_0 */
// Setting: \Blue{$N \gg 1$}, 
// large matrix \Blue{$\VA\in\bbK^{n,n}$}
for(int j = 0; j < N; ++j){
	x = A.lu().solve(b);
	b = some_function(x);
}
/* SAM_LISTING_END_0 */


return 0;
}
