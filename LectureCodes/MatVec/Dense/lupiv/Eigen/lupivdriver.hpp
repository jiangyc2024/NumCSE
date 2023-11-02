///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

namespace lupiv {


using Eigen::MatrixXd;
using Eigen::Upper;
using Eigen::StrictlyLower;

inline
/* SAM_LISTING_BEGIN_0 */
void lupivdriver(const MatrixXd &A, MatrixXd &L, MatrixXd &U){
	MatrixXd Adec = A;
	lupiv(Adec);
	U = Adec.triangularView<Upper>();
	L.setIdentity();
	L += Adec.triangularView<StrictlyLower>();
}
/* SAM_LISTING_END_0 */


} //namespace lupiv