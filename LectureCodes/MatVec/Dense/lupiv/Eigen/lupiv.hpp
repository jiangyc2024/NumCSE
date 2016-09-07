#pragma once

#include <limits>

#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
void lupiv(MatrixXd &A){//insitu
	int n = A.rows();
	for(int k = 0; k < n-1; ++k){
		int j; double p; // p = relativly largest pivot, j = pivot row index 
		p = ( A.col(k).tail(n-k).cwiseAbs().cwiseQuotient( A.block(k,k,n-k,n-k).cwiseAbs().rowwise().maxCoeff() ) ).maxCoeff(&j); // \label{gsrp:1}
		if( p <  std::numeric_limits<double>::epsilon() * A.block(k,k,n-k,n-k).norm() ) // \label{gsrp:2} 
			throw std::logic_error("nearly singular");
		A.row(k).tail(n-k-1).swap(A.row(k+j).tail(n-k-1));// \label{gsrp:3} 
		VectorXd fac = A.col(k).tail(n-k-1) / A(k,k);// \label{gsrp:f}
		A.bottomRightCorner(n-k-1,n-k-1) -= fac * A.row(k).tail(n-k-1);// \label{gsrp:4}  
		A.col(k).tail(n-k-1) = fac;// \label{gsrp:5} 
	}
}
/* SAM_LISTING_END_0 */
