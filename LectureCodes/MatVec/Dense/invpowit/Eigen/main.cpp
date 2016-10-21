///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
#include "invpowit.hpp"

int main(){	
	int n = 10;
	double tol = 1e-6;
/* SAM_LISTING_BEGIN_0 */
MatrixXd A = MatrixXd::Random(n,n);
MatrixXd B = MatrixXd::Random(n,n);
VectorXd ev = invpowit<VectorXd>(A+B, tol);
/* SAM_LISTING_END_0 */

}
