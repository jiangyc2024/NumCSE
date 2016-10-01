/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <stdexcept>
#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! Solving the overdetermined linear system of equations
//! \Blue{$\VA\Vx=\Vb$} by solving \cor{normal equations} \eqref{eq:normeq}
//! The least squares solution is returned by value
VectorXd normeqsolve(const MatrixXd &A,const VectorXd &b) {
  if (b.size() != A.rows()) throw runtime_error("Dimension mismatch");
  // Cholesky solver
  VectorXd x = (A.transpose()*A).llt().solve(A.transpose()*b);
  return x;
}
/* SAM_LISTING_END_0 */
