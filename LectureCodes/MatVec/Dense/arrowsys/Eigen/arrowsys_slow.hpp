///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

namespace arrowsys {


using Eigen::VectorXd;
using Eigen::MatrixXd;

//NOLINTBEGIN(bugprone-easily-swappable-parameters)
inline
/* SAM_LISTING_BEGIN_0 */
VectorXd arrowsys_slow(const VectorXd &d, 
                       const VectorXd &c,
                       const VectorXd &b, double alpha,
                       const VectorXd &y) {
  const Eigen::Index n = d.size();
  MatrixXd A(n + 1, n + 1); // Empty dense matrix
  A.setZero();              // Initialize with all zeros.
  A.diagonal().head(n) = d; // Initialize matrix diagonal from a vector.
  A.col(n).head(n) = c;     // Set rightmost column $\cob{\Vc}$.
  A.row(n).head(n) = b;     // Set bottom row $\cob{\Vb^{\top}}$.
  A(n, n) = alpha;          // Set bottom-right entry $\cob{\alpha}$.
  return A.lu().solve(y);   // Gaussian elimination
}
/* SAM_LISTING_END_0 */
//NOLINTEND(bugprone-easily-swappable-parameters)


} // namespace arrowsys