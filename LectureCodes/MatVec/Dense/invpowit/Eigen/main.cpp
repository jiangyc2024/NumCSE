///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
#include "invpowit.hpp"

int main(){	
	const int n = 10;
	const double tol = 1e-6;
/* SAM_LISTING_BEGIN_0 */
const MatrixXd A = MatrixXd::Random(n,n);
const MatrixXd B = MatrixXd::Random(n,n);
const auto ev = invpowit<VectorXd>(A+B, tol);
/* SAM_LISTING_END_0 */

}
