#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <figure/figure.hpp>
using namespace std;
using namespace Eigen;

int main () {
	#pragma begin<0>
	MatrixXd A(11,11); A.setIdentity();
	A.col(0).setOnes(); A.row(0) = RowVectorXd::LinSpaced(11,11,1);
	// Permutation matrix ($\to$ Def.~\ref{def:permmat}) encoding cyclic permutation
	MatrixXd P(11,11); P.setIdentity(); P(0,0) = 0; P(10,0) = 1;
	mgl::Figure fig1, fig2, fig3, fig4;
	fig1.spy(A);	fig1.setFontSize(4);
	fig1.grid(false);	fig1.save("InvArrowSpy_cpp");
	fig2.spy((P*A*P.transpose()).eval()); fig2.setFontSize(4);
	fig2.grid(false);	fig2.save("ArrowSpy_cpp");
	#pragma end<0>
	return 0;
}
