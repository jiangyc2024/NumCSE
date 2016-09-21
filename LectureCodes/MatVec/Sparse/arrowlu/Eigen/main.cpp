///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace std;
using namespace Eigen;

int main () {
/* SAM_LISTING_BEGIN_0 */
// Build matrix
MatrixXd A(11,11); A.setIdentity();
A.col(10).setOnes(); A.row(10).setOnes();
// A.reverseInPlace(); // used in\cref{ex:arrowmatrixlu}
auto solver = A.lu();
MatrixXd L = MatrixXd::Identity(11,11);
L += solver.matrixLU().triangularView<StrictlyLower>();
MatrixXd U = solver.matrixLU().triangularView<Upper>();
MatrixXd Ainv = A.inverse();
// Plotting
mgl::Figure fig1, fig2, fig3, fig4;
fig1.spy(A);	fig1.setFontSize(4);
fig1.title("Pattern of A");	fig1.save("Apat_cpp");
fig2.spy(L);	fig2.setFontSize(4);
fig2.title("Pattern of L");	fig2.save("Lpat_cpp");
fig3.spy(U);	fig3.setFontSize(4);
fig3.title("Pattern of U");	fig3.save("Upat_cpp");
fig4.spy(Ainv);	fig4.setFontSize(4);
fig4.title("Pattern of A^{-1}");	fig4.save("Ainvpat_cpp");
/* SAM_LISTING_END_0 */
	return 0;
}
