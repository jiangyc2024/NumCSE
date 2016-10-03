///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! plane (2D) Givens rotation avoiding cancellation
void planerot(const VectorXd& a, MatrixXd& G, VectorXd& x){
	if(a(1) != 0){
		double r = a.norm();
		G.row(0) = a.transpose() / r;
		G(1,0) = - a(1) / r;
		G(1,1) = a(0) / r;
		x(0) = r;	x(1) = 0;
	}
	else
		G.setIdentity();
}
/* SAM_LISTING_END_0 */
