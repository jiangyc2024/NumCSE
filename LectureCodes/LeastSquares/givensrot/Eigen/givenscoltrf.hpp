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

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
// Orthogonal transformation of a (column) vector into a multiple of
// the first unit vector by successive Givens transformations
// Note that the output vector could be computed much more efficiently!
void givenscoltrf(const VectorXd &aIn, MatrixXd &Q, VectorXd &aOut) {
  const int n = aIn.size();
  // Assemble rotations in a dense matrix $\cob{\VQ}$
  // For (more efficient) alternatives see Rem.~\cref{rem:storeQ}
  Q.setIdentity(); // Start from $\cob{\VQ=\VI}$
  Matrix2d G;
  aOut = aIn;
  for (int j = 1; j < n; ++j) {
    double a0 = aOut[0], a1 = aOut[j];
    // Determine entries of 2D rotation matrix, see Code~\ref{cpp:planerot}
    double s, c; // s $\leftrightarrow$ $\cob{\sigma}$, c $\leftrightarrow$ $\cob{\gamma}$
    if (a1 != 0.0) {
      if (std::abs(a1) > std::abs(a0)) { // Avoid cancellation/overflow
        double t = -a0 / a1;
        s = 1.0 / std::sqrt(1.0 + t * t);
        c = s * t;
      } else {
        double t = -a1 / a0;
        c = 1.0 / std::sqrt(1.0 + t * t);
        s = c * t;
      }
      G << c, s, -s, c; // Form $2\times 2$ Givens rotation matrix
    } else {            // No rotation required
      G.setIdentity();
    }
    // select 1st and jth element of aOut and use the Map function
    // to prevent copying; equivalent to aOut([1,j]) in \matlab
    Map<VectorXd, 0, InnerStride<>> aOutMap(aOut.data(), 2, InnerStride<>(j));
    aOutMap = G.transpose() * aOutMap;
    // select 1st and jth column of Q (Q(:,[1,j]) in \matlab)
    Map<MatrixXd, 0, OuterStride<>> QMap(Q.data(), n, 2, OuterStride<>(j * n));
    // Accumulate orthogonal transformations in a dense matrix; just done for
    // demonstration purposes! See \cref{rem:storeQ}
    QMap = QMap * G;
  }
}
/* SAM_LISTING_END_0 */
