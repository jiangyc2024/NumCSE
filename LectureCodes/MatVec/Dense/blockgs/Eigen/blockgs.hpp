///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

namespace blockgs {
  

inline
/* SAM_LISTING_BEGIN_0 */
//! in-situ Gaussian elimination, no pivoting
//! right hand side in rightmost column of \Blue{$\VA$}
//! back substitution is not done in this code!
void blockgs(Eigen::MatrixXd &A){
  const Eigen::Index n = A.rows();
  for(Eigen::Index i = 1; i < n; ++i){
    // \Red{rank-1 modification} of \Blue{$\VC$}
    A.bottomRightCorner(n-i,n-i+1) -= A.col(i-1).tail(n-i) * A.row(i-1).tail(n-i+1) / A(i-1,i-1);
    A.col(i - 1).tail(n - i).setZero(); // set $\Vd=0$ \Label[line]{bgscpp:1}
  }
}
/* SAM_LISTING_END_0 */


} //namespace blockgs