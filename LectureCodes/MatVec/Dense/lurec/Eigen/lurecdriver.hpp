///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! post-processing: extract \Blue{$\VL$} and \Blue{$\VU$}
void lurecdriver(const MatrixXd &A, MatrixXd &L, MatrixXd &U){
	MatrixXd A_dec = lurec(A);
	// post-processing: 
	//extract \Blue{$\VL$} and \Blue{$\VU$}
	U = A_dec.triangularView<Upper>();
	L.setIdentity();
	L += A_dec.triangularView<StrictlyLower>();
}
/* SAM_LISTING_END_0 */
