///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>
#include <cassert>
#include <numeric>

/* SAM_LISTING_BEGIN_0 */
//! Computation of $\Vy = \operatorname{triu}(\VA\VB^{T})\Vx$
//! Efficient implementation with backward cumulative sum
//! (partial\_sum)
template<class Vec, class Mat>
void lrtrimulteff(const Mat& A, const Mat& B, const Vec& x, Vec& y){
  const int n = A.rows();
  const int p = A.cols();
  assert( n == B.rows() && p == B.cols()); // size mismatch
  for(int l = 0; l < p; ++l){
    Vec tmp = (B.col(l).array() * x.array()).matrix().reverse();
    std::partial_sum(tmp.begin(), tmp.end(), tmp.begin());
    y += (A.col(l).array() * tmp.reverse().array()).matrix();
  }
}
/* SAM_LISTING_END_0 */
