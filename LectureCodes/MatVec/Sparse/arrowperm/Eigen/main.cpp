#include <iostream>

#include <Eigen/Dense>
#include <figure/figure.hpp>

using namespace std;
using namespace Eigen;

int main () {
/* SAM_LISTING_BEGIN_0 */
MatrixXd A(11,11); A.setIdentity();
A.col(0).setOnes(); A.row(0) = RowVectorXd::LinSpaced(11,11,1);
// Permutation matrix ($\to$ Def.~\ref{def:permmat}) encoding cyclic permutation
MatrixXd P(11,11); P.setZero();
P.topRightCorner(10,10).setIdentity(); P(10,0) = 1;
mgl::Figure fig1, fig2;
fig1.spy(A);	fig1.setFontSize(4);
fig1.save("InvArrowSpy_cpp");
fig2.spy((P*A*P.transpose()).eval()); fig2.setFontSize(4);
fig2.save("ArrowSpy_cpp");
/* SAM_LISTING_END_0 */
	return 0;
}
