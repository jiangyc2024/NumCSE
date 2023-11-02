///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

namespace lurec {


using Eigen::MatrixXd;
using Eigen::VectorXd;

inline
/* SAM_LISTING_BEGIN_0 */
//! in situ recursive LU-factorization
MatrixXd lurec(const MatrixXd &A){
	const Eigen::Index n = A.rows();
	MatrixXd result(n,n);
	if(n > 1){
		const VectorXd fac = A.col(0).tail(n-1) / A(0,0);//\label{lurec:1}
		result.bottomRightCorner(n-1,n-1) = lurec( A.bottomRightCorner(n-1,n-1)	- fac * A.row(0).tail(n-1) );//\label{lurec:2}
		result.row(0) = A.row(0); result.col(0).tail(n-1) = fac;
		return result;
	}
	return A;
}
/* SAM_LISTING_END_0 */


} //namespace lurec