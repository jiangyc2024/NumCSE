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

namespace smw {


using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Lower;
using Eigen::Upper;

/* SAM_LISTING_BEGIN_0 */
// Solving rank-1 updated LSE based on \eqref{eq:smwx}
template <class LUDec>
Eigen::VectorXd smw( const Eigen::VectorXd &u, const Eigen::VectorXd &v, 
    const LUDec &lu, const Eigen::VectorXd &b) {
  const double singfac = 1.0E6; // Do not lose more than 10 digits 
  const Eigen::VectorXd z = lu.solve(b); // $\cob{\Vz=\VA^{-1}\Vb}$ \Label[line]{smw:1}
  const Eigen::VectorXd w = lu.solve(u); // $\cob{\Vw=\VA^{-1}\Vu}$ \Label[line]{smw:2}
  const double alpha = 1.0 + v.dot(w); // Compute denominator of \eqref{eq:smwx}
  const double beta = v.dot(z);        // Factor for numerator of \eqref{eq:smwx}
  if (std::abs(alpha) <= singfac * std::numeric_limits<double>::epsilon()) {
    throw std::runtime_error("A nearly singular");
  }
  return (z - w * beta / alpha); // see \eqref{eq:smwx}
}
/* SAM_LISTING_END_0 */

// Old version
inline
/* SAM_LISTING_BEGIN_1 */
VectorXd smw_old(const MatrixXd &L, const MatrixXd &U, const VectorXd &u,
                 const VectorXd &v, const VectorXd &b) {
  const VectorXd z = U.triangularView<Upper>().solve(
      L.triangularView<Lower>().solve(b)); //\label{smw:1}
  const VectorXd w = U.triangularView<Upper>().solve(
      L.triangularView<Lower>().solve(u)); //\label{smw:2}
  const double alpha = 1.0 + v.dot(w);
  assert(std::abs(alpha) >
             std::numeric_limits<double>::epsilon() * U.lpNorm<1>() &&
         "A nearly singular");
  return z - w * v.dot(z) / alpha;
}
/* SAM_LISTING_END_1 */


} //namespace smw
