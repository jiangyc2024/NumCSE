#pragma once

#include <iostream>
#include <cmath>
#include <limits>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! Solving an LSE \Blue{$\VA\Vx=\Vb$} by Gaussian elimination with partial pivoting
//! \texttt{A} must be an \Blue{$n\times n$}-matrix, \texttt{b} an \Blue{$n$}-vector
void gepiv(const MatrixXd &A, const VectorXd& b, VectorXd& x){
	int n = A.rows();
	MatrixXd Ab(n,n+1);
	Ab << A, b;	//\label{lup:0}
	// Forward elimination by rank-1 modification, see Rem.~\ref{rem:blockgs}
	for(int k = 0; k < n-1; ++k){
		int j; double p; // p = relativly largest pivot, j = pivot row index 
		p = ( Ab.col(k).tail(n-k).cwiseAbs().cwiseQuotient(	Ab.block(k,k,n-k,n-k).cwiseAbs().rowwise().maxCoeff() ) ).maxCoeff(&j);// \label{lup:1}
		if( p <  std::numeric_limits<double>::epsilon() * Ab.block(k,k,n-k,n-k).norm() )
			throw std::logic_error("nearly singular");// \label{lup:2} 
		Ab.row(k).tail(n-k+1).swap(Ab.row(k+j).tail(n-k+1));// \label{lup:3} 
		Ab.bottomRightCorner(n-k-1,n-k) -= Ab.col(k).tail(n-k-1) * Ab.row(k).tail(n-k) / Ab(k,k);// \label{lup:4} 
	}
	// \Hyperlink{RUECKSUBST}{Back substitution} (same as in Code~\ref{cpp:gausselimsolve})
	Ab(n-1,n) = Ab(n-1,n) / Ab(n-1,n-1);
	for(int i = n-2; i >= 0; --i){
		for(int l = i+1; l < n; ++l){
			Ab(i,n) -= Ab(l,n)*Ab(i,l);
		}
		Ab(i,n) /= Ab(i,i);
	}
	x = Ab.rightCols(1); // \label{lup:last}
}
/* SAM_LISTING_END_0 */
