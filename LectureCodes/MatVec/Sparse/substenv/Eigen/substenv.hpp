///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

namespace envelope {


using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;

inline
/* SAM_LISTING_BEGIN_0 */
//! evelope aware forward substitution for \Blue{$\VL\Vx=\Vy$}
//! (\Blue{$\VL$} = lower triangular matrix)
//! argument \texttt{mr}: row bandwidth vector
VectorXd substenv(const MatrixXd &L, const VectorXd &y, const VectorXi &mr) {
  const Eigen::Index n = L.cols();
  VectorXd x(n);
  x(0) = y(0) / L(0, 0);
  for (Eigen::Index i = 1; i < n; ++i) {
    if (mr(i) > 0) {
      const double zeta =
          L.row(i).segment(i - mr(i), mr(i)) * x.segment(i - mr(i), mr(i));
      x(i) = (y(i) - zeta) / L(i, i);
    } else {
      x(i) = y(i) / L(i, i);
    }
  }
  return x;
}
/* SAM_LISTING_END_0 */


} //namespace envelope