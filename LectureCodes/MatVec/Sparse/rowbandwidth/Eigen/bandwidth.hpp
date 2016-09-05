#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! computes rowbandwidth numbers \Blue{$m_i^R(\VA)$} of \Blue{$\VA$} (sparse matrix) according to \cref{def:envelope}
template <class numeric_t>
VectorXi rowbandwidth(const SparseMatrix<numeric_t> &A){
	VectorXi m = VectorXi::Zero(A.rows());
	for(int k = 0; k < A.outerSize(); ++k)
		for(typename SparseMatrix<numeric_t>::InnerIterator it(A,k); it; ++it)
			m(it.row()) = std::max<VectorXi::Scalar>(m( it.row() ), it.row()-it.col() );
	return m;
}
//! computes row bandwidth numbers \Blue{$m_i^R(\VA)$} of \Blue{$\VA$} (dense matrix) according to \cref{def:envelope}
template <class Derived>
VectorXi rowbandwidth(const MatrixBase<Derived> &A){
	VectorXi m = VectorXi::Zero(A.rows());
	for(int i = 1; i < A.rows(); ++i)
		for(int j = 0; j < i; ++j)
			if(A(i,j) != 0){
				m(i) = i - j;
				break;
			}
	return m;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
//! computes colbandwidth numbers \Blue{$m_i^R(\VA)$} of \Blue{$\VA$} (sparse matrix) according to \cref{def:envelope}
template <class numeric_t>
VectorXi colbandwidth(const SparseMatrix<numeric_t> &A){
	VectorXi m = VectorXi::Zero(A.cols());
	for(int k = 0; k < A.outerSize(); ++k)
		for(typename SparseMatrix<numeric_t>::InnerIterator it(A,k); it; ++it)
			m(it.col()) = std::max<VectorXi::Scalar>(m( it.col() ), it.col()-it.row() );
	return m;
}


//! computes column bandwidth numbers \Blue{$m_i^R(\VA)$} of \Blue{$\VA$} (dense matrix) according to \cref{def:envelope}
template <class Derived>
VectorXi colbandwidth(const MatrixBase<Derived> &A){
	VectorXi m = VectorXi::Zero(A.cols());
	for(int j = 1; j < A.cols(); ++j)
		for(int i = 0; i < j; ++i)
			if(A(i,j) != 0){
				m(j) = j - i;
				break;
			}
	return m;
}
/* SAM_LISTING_END_1 */
