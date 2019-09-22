///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
///            R. Hiptmair <hiptmair@sam.math.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cassert>
#include <iostream>
#include <limits>

#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
// Solving rank-1 updated LSE based on \eqref{eq:smwx}
template <class LUDec>
VectorXd smw(const LUDec &lu, const VectorXd &u, const VectorXd &v,
             const VectorXd &b) {
  const VectorXd z = lu.solve(b); // $\cob{\Vz=\VA^{-1}\Vb}$ \Label[line]{smw:1}
  const VectorXd w = lu.solve(u); // $\cob{\Vw=\VA^{-1}\Vu}$ \Label[line]{smw:2}
  double alpha = 1.0 + v.dot(w);  // Compute denominator of \eqref{eq:smwx}
  double vdz = v.dot(z);          // Factor for numerator of \eqref{eq:smwx}
  if (std::abs(alpha) < std::numeric_limits<double>::epsilon() * std::abs(vdz))
    throw std::runtime_error("A nearly singular");
  else
    return (z - w * vdz / alpha); // see \eqref{eq:smwx}
}
/* SAM_LISTING_END_0 */

// Old version
/* SAM_LISTING_BEGIN_1 */
VectorXd smw_old(const MatrixXd &L, const MatrixXd &U, const VectorXd &u,
                 const VectorXd &v, const VectorXd &b) {
  VectorXd z = U.triangularView<Upper>().solve(
      L.triangularView<Lower>().solve(b)); //\label{smw:1}
  VectorXd w = U.triangularView<Upper>().solve(
      L.triangularView<Lower>().solve(u)); //\label{smw:2}
  double alpha = 1.0 + v.dot(w);
  assert(std::abs(alpha) >
             std::numeric_limits<double>::epsilon() * U.lpNorm<1>() &&
         "A nearly singular");
  return z - w * v.dot(z) / alpha;
}
/* SAM_LISTING_END_1 */
