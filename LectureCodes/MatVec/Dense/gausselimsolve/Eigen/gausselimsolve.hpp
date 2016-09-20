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
//! Gauss elimination without pivoting, \texttt{x = A\symbol{92}b}
//! \texttt{A} must be an \Blue{$n\times n$}-matrix, \texttt{b} an \Blue{$n$}-vector
void gausselimsolve(const MatrixXd &A, const VectorXd& b, VectorXd& x){
	int n = A.rows();
	MatrixXd Ab(n,n+1);
	Ab << A, b;//\label{cppgse:1}
	// Forward elimination (\textit{cf.} step \ding{192} in Ex.~\ref{ex:GE})
	for(int i = 0; i < n-1; ++i){
		double pivot = Ab(i,i);
		for(int k = i+1; k < n; ++k){
			double fac = Ab(k,i)/pivot;
			Ab.block(k,i+1,1,n-i)-= fac * Ab.block(i,i+1,1,n-i);//\label{cppgse:vec}
		}
	}
	// \Hyperlink{RUECKSUBST}{Back substitution} (\textit{cf.} step \ding{193} in Ex.~\ref{ex:GE})
	Ab(n-1,n) = Ab(n-1,n) / Ab(n-1,n-1);
	for(int i = n-2; i >= 0; --i){
		for(int l = i+1; l < n; ++l){
			Ab(i,n) -= Ab(l,n)*Ab(i,l);
		}
		Ab(i,n) /= Ab(i,i);
	}
	x = Ab.rightCols(1);// \label{cppgse:last}
}
/* SAM_LISTING_END_0 */
