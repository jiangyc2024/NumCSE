#pragma once

#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
void lupivdriver(const MatrixXd &A, MatrixXd &L, MatrixXd &U){
	MatrixXd Adec = A;
	lupiv(Adec);
	U = Adec.triangularView<Upper>();
	L.setIdentity();
	L += Adec.triangularView<StrictlyLower>();
}
/* SAM_LISTING_END_0 */
