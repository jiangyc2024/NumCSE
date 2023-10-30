///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Dense>

namespace gepiv {


using Eigen::MatrixXd;
using Eigen::VectorXd;

inline
/* SAM_LISTING_BEGIN_0 */
//! Solving an LSE \Blue{$\VA\Vx=\Vb$} by Gaussian elimination with partial pivoting
//! \texttt{A} must be an \Blue{$n\times n$}-matrix, \texttt{b} an \Blue{$n$}-vector
void gepiv(const MatrixXd &A, const VectorXd& b, VectorXd& x){
	const Eigen::Index n = A.rows();
	MatrixXd Ab(n,n+1);
	Ab << A, b;	//\label{lup:0}
	// Forward elimination by rank-1 modification, see Rem.~\ref{rem:blockgs}
	for(Eigen::Index k = 0; k < n-1; ++k){
		Eigen::Index j = -1; // j = pivot row index 
		// p = relatively largest pivot
		const double p = ( Ab.col(k).tail(n-k).cwiseAbs().cwiseQuotient(	Ab.block(k,k,n-k,n-k).cwiseAbs().rowwise().maxCoeff() ) ).maxCoeff(&j);// \label{lup:1}
		if( p <  std::numeric_limits<double>::epsilon() * Ab.block(k,k,n-k,n-k).norm() ) {
			throw std::logic_error("nearly singular");// \label{lup:2} 
		}
		Ab.row(k).tail(n-k+1).swap(Ab.row(k+j).tail(n-k+1));// \label{lup:3} 
		Ab.bottomRightCorner(n-k-1,n-k) -= Ab.col(k).tail(n-k-1) * Ab.row(k).tail(n-k) / Ab(k,k);// \label{lup:4} 
	}
	// \Hyperlink{RUECKSUBST}{Back substitution} (same as in Code~\ref{cpp:gausselimsolve})
	Ab(n-1,n) = Ab(n-1,n) / Ab(n-1,n-1);
	for(Eigen::Index i = n-2; i >= 0; --i){
		for(Eigen::Index l = i+1; l < n; ++l){
			Ab(i,n) -= Ab(l,n)*Ab(i,l);
		}
		Ab(i,n) /= Ab(i,i);
	}
	x = Ab.rightCols(1); // \label{lup:last}
}
/* SAM_LISTING_END_0 */


} //namespace gepiv
