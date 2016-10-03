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

#include "planerot.hpp"

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! Orthogonal transformation of a (column) vector into a multiple of the first unit vector by successive Givens transformations
void givenscoltrf(const VectorXd& aIn, MatrixXd& Q, VectorXd& aOut){
	unsigned int n = aIn.size();
	Q.setIdentity(); // Assemble rotations in matrix, alternative see Rem.~\ref{rem:orthstore}
	MatrixXd G(2,2);
	VectorXd tmp(2), xDummy(2);
	aOut = aIn;
	for(int j = 1; j < n; ++j){
		tmp(0) = aOut(0); tmp(1) = aOut(j);
		planerot(tmp, G, xDummy); // see Code~\ref{cpp:planerot}
		Map<VectorXd, 0, InnerStride<> > aOutMap(aOut.data(), 2, InnerStride<>(j)); // select 1st and jth element of aOut and use the Map function to prevent copying, equivalent to aOut([1,j]) in \matlab
		aOutMap = G * aOutMap;
		Map<MatrixXd, 0, OuterStride<> > QMap(Q.data(), n, 2, OuterStride<>(j*n)); // select 1st and jth column of Q, equivalent to Q(:,[1,j]) in \matlab
		QMap = QMap * G.transpose();
	}
}
/* SAM_LISTING_END_0 */
