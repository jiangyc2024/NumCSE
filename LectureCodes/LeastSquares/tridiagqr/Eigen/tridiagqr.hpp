///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>
#include <Eigen/Dense>

#include "planerot.hpp"

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! @brief Solves the tridiagonal system \Blue{$\VA\Vx=\Vb$} with QR-decomposition
//! @param[in] \Blue{$\Vd$} Vector of dim $n$; the diagonal elements
//! @param[in] \Blue{$\Vc$} Vector of dim $n-1$; the lower diagonal elements
//! @param[in] \Blue{$\Ve$} Vector of dim $n-1$; the upper diagonal elements
//! @param[in] \Blue{$\Vb$} Vector of dim $n$; the rhs.
//! @param[out] \Blue{$\Vx$} Vector of dim $n$
VectorXd tridiagqr(VectorXd c, VectorXd d, VectorXd e, VectorXd& b){
  int n = d.size();
  // resize the vectors c and d to correct length if needed
  c.conservativeResize(n); e.conservativeResize(n);
  double t = d.norm() + e.norm() + c.norm();
  Matrix2d R;	Vector2d z, tmp;
  for(int k = 0; k < n-1; ++k){
    tmp(0) = d(k); tmp(1) = e(k);
    // Use givensrotation to set the entries below the diagonal
    // to zero
    planerot(tmp, R, z); // see Code~\ref{cpp:planerot}
    if( std::abs(z(0))/t < std::numeric_limits<double>::epsilon() )
      throw std::runtime_error("A nearly singular");
    // Update all other entries of the matrix and rhs. which 
    // were affected by the givensrotation
    d(k) = z(0);
    b.segment(k,2).applyOnTheLeft(R);	// rhs.
    // Block of the matrix affected by the givensrotation
    Matrix2d Z;
    Z << c(k), 	0,  d(k+1), c(k+1);
    Z.applyOnTheLeft(R);
    // Write the transformed block back to the corresponding places
    c.segment(k,2) = Z.diagonal(); d(k+1) = Z(1,0); e(k) = Z(0,1);
  }
  // Note that the \Blue{$\Ve$} is now above \Blue{$\Vd$} and \Blue{$\Vc$}
  // Backsubstitution acting on upper triangular matrix
  // with upper bandwidth 2 (stored in vectors).
  VectorXd x(n);
  // last row
  x(n-1) = b(n-1)/d(n-1);
  if(n >= 2) {
    // 2nd last row
    x(n-2) = (b(n-2)-c(n-2)*x(n-1))/d(n-2);
    // remaining rows
    for(int i = n-3; i >= 0; --i)
      x(i) = ( b(i) - c(i) * x(i+1) - e(i)*x(i+2) ) / d(i);
  }
  return x;
}
/* SAM_LISTING_END_0 */
