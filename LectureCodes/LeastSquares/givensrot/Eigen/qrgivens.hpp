///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>
#include <cmath>

#include "planerot.hpp"

namespace qrgivens {


using Eigen::MatrixXd;
using Eigen::Matrix2d;
using Eigen::Vector2d;
using planerot::planerot;

inline
/* SAM_LISTING_BEGIN_0 */
//! QR decomposition of \emph{square} matrix \texttt{A} by successive Givens
//! transformations
void qrgivens(const MatrixXd &A, MatrixXd &Q, MatrixXd &R) {
  const Eigen::Index n = A.rows();
  // Assemble rotations in a dense matrix.
  // For (more efficient) alternatives see Rem.~\cref{rem:storeQ}
  Q.setIdentity();
  Matrix2d G;
  Vector2d tmp;
  Vector2d xDummy;
  R = A; // In situ transformation
  for (Eigen::Index i = 0; i < n - 1; ++i) {
    for (Eigen::Index j = n - 1; j > i; --j) {
      tmp(0) = R(j - 1, i);
      tmp(1) = R(j, i);
      planerot(tmp, G, xDummy); // see Code~\ref{cpp:planerot}
      R.block(j - 1, 0, 2, n) = G.transpose() * R.block(j - 1, 0, 2, n);
      Q.block(0, j - 1, n, 2) = Q.block(0, j - 1, n, 2) * G;
    }
  }
}
/* SAM_LISTING_END_0 */


} // namespace qrgivens