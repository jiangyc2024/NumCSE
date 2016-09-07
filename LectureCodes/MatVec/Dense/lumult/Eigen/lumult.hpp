#pragma once

#include <cassert>

#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! Multiplication of normalized lower/upper triangular matrices
void lumult(const MatrixXd& L, const MatrixXd& U, MatrixXd& A){
  int n = L.rows();
  assert( n == L.cols() && n == U.cols() && n == U.rows() );
  A = MatrixXd::Zero(n,n);
  for(int k = 0; k < n; ++k){
  	for(int j = k; j < n; ++j)
      A(k,j) = U(k,j) + (L.block(k,0,1,k) * U.block(0,j,k,1))(0,0);
  	for(int i = k+1; i < n; ++i)
  	  A(i,k) = (L.block(i,0,1,k) * U.block(0,k,k,1))(0,0) + L(i,k)*U(k,k);
  }
}
/* SAM_LISTING_END_0 */
