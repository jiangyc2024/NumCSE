#pragma once
//! Computation of \Vy = \operatorname{triu}(\VA\VB^{T})\Vx
//! Efficient implementation with backward cumulative sum(partial_sum)
template<class Vec, class Mat>
void lrtrimulteff(const Mat& A, const Mat& B, const Vec& x, Vec& y){
  int n = A.rows(), p = A.cols();
  assert( n == B.rows() && p == B.cols()); // size missmatch
  for(int l = 0; l < p; ++l){
    Vec tmp = (B.col(l).array() * x.array()).matrix().reverse();
    std::partial_sum(tmp.data(), tmp.data()+n, tmp.data());
    y += (A.col(l).array() * tmp.reverse().array()).matrix();
  }
}
