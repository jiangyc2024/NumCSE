///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
VectorXd arrowsys_slow(const VectorXd &d, const VectorXd &c, const VectorXd &b,
                       const double alpha, const VectorXd &y) {
  int n = d.size();
  MatrixXd A(n + 1, n + 1); // Empty dense matrix
  A.setZero();              // Initialize with all zeros.
  A.diagonal().head(n) = d; // Initializee matrix diagonal from a vector.
  A.col(n).head(n) = c;     // Set rightmost column $\cob{\Vc}$.
  A.row(n).head(n) = b;     // Set bottom row $\cob{\Vb^{\top}}$.
  A(n, n) = alpha;          // Set bottom-right entry $\cob{\alpha}$.
  return A.lu().solve(y);   // Gaussian elimination
}
/* SAM_LISTING_END_0 */
