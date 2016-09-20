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
//! Gauss elimination without pivoting, $\VX = \VA^{-1}\VB$
//! \texttt{A} must be an \Blue{$n\times n$}-matrix, \texttt{B} an \Blue{$n\times m$}-matrix
void gausselimsolvemult(const MatrixXd &A, const MatrixXd& B, MatrixXd& X){
	int n = A.rows();
	int m = B.cols();
	MatrixXd AB(n, n+m);
	AB << A, B;
	// Forward elimination, do not forget the B part of the Matrix
	for(int i = 0; i < n-1; ++i){
		double pivot = AB(i,i);
		for(int k = i+1; k < n; ++k){
			double fac = AB(k,i)/pivot;
			AB.block(k,i+1,1,m+n-i-1)-= fac * AB.block(i,i+1,1,m+n-i-1);
		}
	}
	// Back substitution
	AB.block(n-1, n,1, m) /= AB(n-1,n-1);
	for(int i = n-2; i >= 0; --i){
		for(int l = i+1; l < n; ++l){
			AB.block(i,n,1,m) -= AB.block(l,n,1,m)*AB(i,l);
		}
		AB.block(i,n,1,m) /= AB(i,i);
	}
	X = AB.rightCols(m);
}
/* SAM_LISTING_END_0 */
